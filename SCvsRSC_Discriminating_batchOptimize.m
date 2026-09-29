%% SCvsRSC_Discriminating_batchOptimize.m
% 목적:
% 1) 30분 LFP에서 random-start 30초 sample 10개 = batch 1개
% 2) batch 10개 = 총 100 sample
% 3) 100개 sample에서 boundaryDepth_oneSwitch가 가장 단일값으로 수렴하도록
%    RSC score / SC score weights를 업데이트

clear; clc; close all;

binPath  = "260415_3_g1_tcat.imec0.lf.bin";
metaPath = "260415_3_g1_tcat.imec0.lf.meta"; %#ok<NASGU>

%% 0. 전체 기록 길이 확인
[dataPath, nameOnly, extOnly] = fileparts(binPath);
if strlength(string(dataPath)) == 0
    dataPath = pwd;
end
binName = string(nameOnly) + string(extOnly);
meta = SGLX_readMeta.ReadMeta(char(binName), char(dataPath));
fullDurationSec = str2double(meta.fileTimeSecs);

if fullDurationSec < 30*60
    warning("meta.fileTimeSecs = %.2f sec입니다. 30분보다 짧습니다. 실제 길이를 기준으로 sampling합니다.", fullDurationSec);
end

%% 1. 30초 random sample 10개 x batch 10개 생성
rng(1);  % 재현성. 매번 새 random sampling을 원하면 rng('shuffle')로 변경

sampleSec = 30;
nSamplePerBatch = 10;
nBatch = 10;
nTotalSample = nSamplePerBatch * nBatch;

maxStartSec = fullDurationSec - sampleSec;
if maxStartSec <= 0
    error("The file is shorter than sampleSec = %.1f sec.", sampleSec);
end

startSecMat = rand(nSamplePerBatch, nBatch) * maxStartSec;
startSecVec = startSecMat(:);
batchIDVec = repelem((1:nBatch).', nSamplePerBatch);

%% 2. 각 30초 sample에서 CSD 계산
% 주의: 이 단계가 가장 오래 걸립니다. 각 sample별 CSD는 weight와 무관하므로 한 번만 계산합니다.
csdOut = cell(nTotalSample, 1);

for i = 1:nTotalSample
    opts = struct();
    opts.startSec = startSecVec(i);
    opts.durationSec = sampleSec;
    opts.makeFigure = false;

    fprintf("\n[CSD] sample %03d/%03d | batch %02d | start = %.2f sec\n", ...
        i, nTotalSample, batchIDVec(i), opts.startSec);

    csdOut{i} = CSD_Calculating_batchable(binPath, metaPath, opts);
end

%% 3. 안정화할 score weight 탐색
% weight 정의:
% RSC = [width, internalCorr, internalGradient, thetaGamma, boundarySharpness]
% SC  = [amplitude, sinkSourceAsymmetry, localPeak, internalGradient]
% depthPrior = scalar

baseWeights.RSC = [1.00, 1.00, -1.00, 0.70, -0.50];
baseWeights.SC  = [1.00, 1.00, 0.70, 0.50];
baseWeights.depthPrior = 1.50;

% 탐색 범위: 기존 weight 주변에서 multiplicative/random perturbation
nCandidate = 50;
boundaryTolUm = 20;  % 같은 boundary로 볼 depth tolerance. Neuropixels spacing에 맞춰 20 또는 40 권장

candidates = repmat(baseWeights, nCandidate + 1, 1);
candidates(1) = baseWeights;

for c = 2:nCandidate+1
    scaleRSC = 0.5 + 1.5 * rand(size(baseWeights.RSC));  % 0.5~2.0배
    scaleSC  = 0.5 + 1.5 * rand(size(baseWeights.SC));

    candidates(c).RSC = baseWeights.RSC .* scaleRSC;
    candidates(c).SC  = baseWeights.SC  .* scaleSC;
    candidates(c).depthPrior = 0.5 + 2.5 * rand();       % 0.5~3.0

    % sign constraint: 기존 해석 방향은 유지
    candidates(c).RSC(3) = -abs(candidates(c).RSC(3)); % RSC: low internalGradient 선호
    candidates(c).RSC(5) = -abs(candidates(c).RSC(5)); % RSC: low boundarySharpness 선호
    candidates(c).SC(4)  = +abs(candidates(c).SC(4));  % SC: high internalGradient 선호
end

optTable = table();
bestObj = inf;
bestWeights = baseWeights;
bestBoundaryVec = nan(nTotalSample, 1);

for c = 1:numel(candidates)
    boundaryVec = nan(nTotalSample, 1);

    for i = 1:nTotalSample
        r = SCvsRSC_single_from_out(csdOut{i}, candidates(c), false, false);
        boundaryVec(i) = r.boundaryDepth_oneSwitch;
    end

    validBoundary = boundaryVec(isfinite(boundaryVec));

    if isempty(validBoundary)
        obj = inf;
        nUnique = inf;
        spreadUm = inf;
        madUm = inf;
        modeBoundary = NaN;
        convergenceRate = 0;
    else
        roundedBoundary = round(validBoundary / boundaryTolUm) * boundaryTolUm;
        u = unique(roundedBoundary);
        counts = arrayfun(@(x) sum(roundedBoundary == x), u);
        [maxCount, maxPos] = max(counts);

        nUnique = numel(u);
        spreadUm = max(validBoundary) - min(validBoundary);
        madUm = mad(validBoundary, 1);
        modeBoundary = u(maxPos);
        convergenceRate = maxCount / numel(validBoundary);

        % 목적함수: unique boundary 수, spread, MAD를 동시에 최소화
        % convergenceRate가 높을수록 감점
        obj = 1000*(nUnique - 1) + 2*spreadUm + 5*madUm - 500*convergenceRate;
    end

    newRow = table(c, obj, nUnique, spreadUm, madUm, modeBoundary, convergenceRate, ...
        'VariableNames', {'candidateID','objective','nUnique','spreadUm','madUm','modeBoundaryUm','convergenceRate'});
    optTable = [optTable; newRow]; %#ok<AGROW>

    if obj < bestObj
        bestObj = obj;
        bestWeights = candidates(c);
        bestBoundaryVec = boundaryVec;
    end

    fprintf("[OPT] %03d/%03d | obj %.2f | nUnique %d | spread %.1f um | mode %.1f um | conv %.2f\n", ...
        c, numel(candidates), obj, nUnique, spreadUm, modeBoundary, convergenceRate);
end

%% 4. 최종 weight로 batch/sample별 boundary 요약
roundedBest = round(bestBoundaryVec / boundaryTolUm) * boundaryTolUm;
finalBoundaryDepth_oneSwitch = mode(roundedBest(isfinite(roundedBest)));

resultTable = table((1:nTotalSample).', batchIDVec, startSecVec, bestBoundaryVec, roundedBest, ...
    'VariableNames', {'sampleID','batchID','startSec','boundaryDepth_oneSwitch','roundedBoundaryDepth'});

batchSummary = groupsummary(resultTable, 'batchID', {'median','std'}, 'boundaryDepth_oneSwitch');

fprintf("\n================ FINAL ================\n");
fprintf("Final converged boundaryDepth_oneSwitch = %.1f um\n", finalBoundaryDepth_oneSwitch);
fprintf("Best RSC weights = [%s]\n", num2str(bestWeights.RSC, ' %.3f'));
fprintf("Best SC  weights = [%s]\n", num2str(bestWeights.SC,  ' %.3f'));
fprintf("Best depthPriorWeight = %.3f\n", bestWeights.depthPrior);

disp(resultTable);
disp(batchSummary);

%% 5. 저장
save("SCvsRSC_batchOptimization_result.mat", ...
    "bestWeights", "bestBoundaryVec", "finalBoundaryDepth_oneSwitch", ...
    "resultTable", "batchSummary", "optTable", "startSecMat", "sampleSec", "boundaryTolUm");

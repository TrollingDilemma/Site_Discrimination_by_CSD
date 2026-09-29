function out = CSD_Calculating_batchable(binPath, metaPath, opts)
%CSD_CALCULATING_BATCHABLE CSD calculation for arbitrary time windows.
% opts.startSec, opts.durationSec, opts.makeFigure

if nargin < 3 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, "startSec"); opts.startSec = 0; end
if ~isfield(opts, "durationSec"); opts.durationSec = []; end
if ~isfield(opts, "makeFigure"); opts.makeFigure = false; end

%% 1. 경로와 파일명
binPath = "C:\Users\DohyunShin\Documents\MATLAB\CSD_demo\260415_3_g1_tcat.imec0.lf.bin"
[dataPath, nameOnly, extOnly] = fileparts(binPath);
if strlength(string(dataPath)) == 0
    dataPath = pwd;
end
binName = string(nameOnly) + string(extOnly);

%% 2. meta 읽기
meta = SGLX_readMeta.ReadMeta(char(binName), char(dataPath));

fs = SGLX_readMeta.SampRate(meta);
nChan = str2double(meta.nSavedChans);

[LF, AP, SY] = SGLX_readMeta.ChannelCountsIM(meta);

fprintf("nSavedChans = %d\n", nChan);
fprintf("Sampling rate = %.2f Hz\n", fs);
fprintf("LF = %d, AP = %d, SY = %d\n", LF, AP, SY);

nNeural = nChan - SY;
chanList = 1:nNeural;

fprintf("Neural LFP channels = %d\n", nNeural);

%% 3. snsGeomMap 파싱

geomText = meta.snsGeomMap;

tokens = regexp(geomText, '\(([^)]*)\)', 'tokens');

geomHeader = tokens{1}{1};
fprintf("snsGeomMap header = %s\n", geomHeader);

nGeom = numel(tokens) - 1;

if nGeom ~= nNeural
    warning("nGeom = %d, nNeural = %d 입니다. geometry와 neural channel 수가 다릅니다.", ...
        nGeom, nNeural);
end

shank = zeros(nGeom, 1);
x_in_um = zeros(nGeom, 1);
y_in_um = zeros(nGeom, 1);
used = zeros(nGeom, 1);

for i = 1:nGeom
    vals = sscanf(tokens{i+1}{1}, '%f:%f:%f:%f');

    shank(i) = vals(1);
    x_in_um(i) = vals(2);
    y_in_um(i) = vals(3);
    used(i) = vals(4);
end

fprintf("Unique shanks: ");
disp(unique(shank)');

fprintf("Unique x positions: ");
disp(unique(x_in_um)');

fprintf("Depth range: %.1f to %.1f um\n", min(y_in_um), max(y_in_um));

%% 4. Gain correction

fullDurationSec = str2double(meta.fileTimeSecs);
if isempty(opts.durationSec)
    durationSec = fullDurationSec;
else
    durationSec = opts.durationSec;
end

if opts.startSec < 0 || opts.startSec + durationSec > fullDurationSec
    error("Requested window [%.3f, %.3f] sec is outside file duration %.3f sec.", ...
        opts.startSec, opts.startSec + durationSec, fullDurationSec);
end

samp0 = round(opts.startSec * fs);
nSamp = round(durationSec * fs);
raw = SGLX_readMeta.ReadBin(samp0, nSamp, meta, char(binName), char(dataPath));
dataV = SGLX_readMeta.GainCorrectIM(raw, chanList, meta);

lfp_V = dataV(chanList, :);
lfp_in_uV = lfp_V * 1e6;

clear raw lfp_V;

%% 5. 시간 방향 preprocessing
lfp_detrend = detrend(double(lfp_in_uV).').'; % detrend는 기본적으로 column 방향으로 작동하므로, transpose -> detrend -> 다시 transpose 형태로 작업

clear lfp_in_uV;

hp = 1;
lp = 300;

[b,a] = butter(3, [hp lp] / (fs/2), "bandpass");
lfp_bp = filtfilt(b, a, lfp_detrend.').'; % filtfilt는 기본적으로 column 방향으로 작동하므로, transpose -> filtering -> 다시 transpose 형태로 작업

clear lfp_detrend;

%% 6. 60Hz line noise filtering

lfp_notch = lfp_bp;
lineFreqs = [60 120 180 240];

for f0 = lineFreqs
    if f0 < fs/2
        wo = f0 / (fs/2);
        bw = wo/ 35;
        [bn, an] = designNotchPeakIIR(CenterFrequency=wo, Bandwidth=bw, Response="notch");
        lfp_notch = filtfilt(bn, an, lfp_notch.').'; % filtfilt는 기본적으로 column 방향으로 작동하므로, transpose -> filtering -> 다시 transpose 형태로 작업
    end
end

clear lfp_bp;

%% 7. Bad channel detection
chanStd = std(lfp_notch, 0 ,2, "omitmissing"); % sample normalization을 row방향으로(즉 시계열로)
medStd = median(chanStd, "omitmissing");

bad = (chanStd > 5*medStd) | (chanStd< 0.2*medStd); % channel recording의 변동성을 기준으로 불량단자 cut-off 
fprintf("Detected bad channels: %d / %d\n", sum(bad), numel(bad));

%% 8. Bad channel interpolation
lfp_interp = lfp_notch;
clear lfp_notch;

allCh  = (1:size(lfp_interp,1)).'; % dimension = [384 X 1 double]
goodCh = find(~bad);
badCh  = find(bad);

if ~isempty(badCh) && numel(goodCh) >= 2
    % 조건1. ~isempty(badCh): interpolation할 대상이 1개 이상 존재한다.
    % 조건2. numel(goodCh) >= 2: interpolation하려면 기준점이 2개 이상은 있어야 한다.
    % 두 조건이 동시에 만족되지 않으면 end
    lfp_interp(badCh,:) = interp1(goodCh, lfp_interp(goodCh,:), ...
        badCh, "linear", "extrap");
end

%% 9. Common median reference
commonRef = median(lfp_interp, 1, "omitnan");
lfp_car = lfp_interp - commonRef;

clear lfp_interp commonRef;

%% 10. meta의 구조 기반으로 depth X time 행렬 생성
targetShank = 0;

validGeomIdx = find(shank == targetShank & used == 1);

depth_in_um = unique(y_in_um(validGeomIdx));
depth_in_um = sort(depth_in_um);

V_depth = nan(numel(depth_in_um), size(lfp_car,2));

for d = 1:numel(depth_in_um)
    thisDepth = depth_in_um(d);

    chAtDepth = validGeomIdx(y_in_um(validGeomIdx) == thisDepth);
    V_depth(d,:) = mean(lfp_car(chAtDepth,:), 1, "omitmissing"); % 같은 depth는 평균치기
end

clear lfp_car;

%% 11. Depth-wise frequency-band feature extraction
% 이로부터 무엇을 할 수 있는가?
% 1) 1500um 근방에서 drastic transition을 detection
% 2) 이후 heatmap에서 분석의 대상에서 제외할 '상층부'를 masking할 수 있음
% 3) RSC-like | SC-like | ambiguous 판정

expectedBoundaryUm = 1500;
boundaryWindowUm = 400;
changeZThreshold = 2.5; % 차후 정의된 boundary의 robust-Z를 계산하였을 때, 이의 마땅함을 평가할 기준 -> 이후 reference 따라 정당한 값으로 수정 가능
excludeSide = "shallower"

% size(V_depth, 1) = (겹치지 않도록 sorting된 depth의 개수)
% size(V_depth, 2) = (durationSec로부터 정의된 분석대상 sample의 개수)
% step for PowerSpectralDensity calculation
winSec = 2;
winLength = round(winSec * fs); % sample#-wise로 length를 정의
winLength = min(winLength, size(V_depth, 2)); % window 길이가 전체 분석 대상보다 길 때를 대비하였음

noverlap = round(0.5 * winLength); % Welch window가 (segment들끼리) 50% 겹치도록 설정
nfft = max(2^10, 2^nextpow2(winLength)); % nfft는 적어도 2^10보다는 크도록 강제

freqBands = struct();
freqBands.delta = [1 4];
freqBands.theta = [4 8];
freqBands.alpha_beta = [8 30];
freqBands.lowGamma = [30 80];
freqBands.highGamma = [80 150];

bandNames = fieldnames(freqBands);
nBand = numel(bandNames);

nDepth = size(V_depth, 1);

meanBandPSD = nan(nDepth, nBand); % uV^2/Hz
bandPower = nan(nDepth, nBand); % uV^2

for d = 1:nDepth
    x = V_depth(d,:);
    x = double(x);
    x = x - mean(x, "omitmissing");

    [Pxx, F] = pwelch(x, winLength, noverlap, nfft, fs);

    for b = 1:nBand
        band = freqBands.(bandNames{b});
        idx = F >= band(1) & F <= band(2);

        if any(idx)
            meanBandPSD(d,b) = mean(Pxx(idx), "omitmissing"); % d번째 depth에서 b번째 frequency band의 평균 PSD
            bandPower(d,b) = trapz(F(idx), Pxx(idx)); % d번째 depth에서 b번째 frequency band의 총 power
        end
    end
end

% dB로 환산
logMeanBandPSD = 10 * log10(meanBandPSD + eps);
logBandPower = 10 * log10(bandPower + eps);

% depth별 RMS 및 adjacent-depth correlation
depthRMS = rms(V_depth, 2, "omitmissing");
logDepthRms = log10(depthRMS + eps);

adjCorr = nan(nDepth - 1, 1);
for d = 1:nDepth-1
    x1 = V_depth(d,:);
    x2 = V_depth(d+1,:);
    
    valid = isfinite(x1) & isfinite(x2);

    if sum(valid) > 10
        x = x1(valid);
        y = x2(valid);
        adjCorr(d) = corr(x(:), y(:), "Rows", "complete");
    else
        adjCorr(d) = NaN;
    end
    % adjCorr(d) = 1: 이웃 depth의 LFP가 거의 같은 시간 패턴으로 움직임
    % adjCorr(d) 낮음: 이웃 depth 사이에서 LFP 파형이 갑자기 달라짐
end
   
corrDrop = 1 - adjCorr; % 인접 depth 사이의 괴리를 설명하는 값

% 변화점 detection용 feature matrix
% feature란? -> RMS + 각 band mean PSD
% 경우1. RMS가 변하지만, correlation이 높을 때-> 구조적 boundary라기보다는, 단순 amplitude gradient일 수 있음
% 경우2. RMS는 일정하지만, correlation이 낮을 때-> 인전 depth 사이에서 신호 source, tissue structure, recording regimen 등이 달라졌음을 시사
featureMat = [logDepthRms, logMeanBandPSD];

% robust z-score normalization: featureMat의 column 각각을 비교 가능할 수 있도록 normalization
% cf) robust Z = (X-median)/MAD -> extreme value의 효과 상쇄
featureMed = median(featureMat, 1, "omitmissing");
featureMAD = mad(featureMat, 1, 1); % y = mad(X,flag,dim)
featureMAD(featureMAD == 0) = eps;

featureZ = (featureMat - featureMed) ./ featureMAD;

% 인접 depth 사이 feature 변화량
featureGradient = sqrt(sum(diff(featureZ,1,1).^2, 2));

% correlation drop을 featureGradient와 합산하여 depth에 따른 transition score로 변환
corrDropZ = (corrDrop - median(corrDrop, "omitmissing")) ./ (mad(corrDrop, 1) + eps);
transitionScore = featureGradient + max(corrDropZ, 0); % max(corrDropZ, 0): corrDropZ와 0행렬을 element-wise로 비교해 큰 값을 반환 -> 양의 evidence만 취함
%% transitionScore 계산 방식에 관하여...
% size(featureGradient) = [191 X 1], size(corrDropZ) = [191 X 1]인 상황에서 transitionScore는 두 column-vector를 element-wise로 더하는 꼴이 된다.
% 이때 커다란 RMS(in Z scale) 자체가 boundary의 가능성을 시사하는 척도이므로, 굳이 negative Z로 감산할 필요는 없다.
% 따라서 max(CorrDropZ, 0)으로 선처리하였음.
%%
midDepth = (depth_in_um(1:end-1) + depth_in_um(2:end)) / 2;

% 1500um 주변에서 boundary 후보 탐색
searchIdx = midDepth >= expectedBoundaryUm - boundaryWindowUm & midDepth <= expectedBoundaryUm + boundaryWindowUm;
% searchIdx의 자료형은 logical, 전체 깊이에 대해 경계 찾는 작업을 진행하고
% 싶으면 searchIdx를

if any(searchIdx)
    localScores = transitionScore(searchIdx);
    localDepths = midDepth(searchIdx);
    [maxLocalScore, maxIdx] = max(localScores);
    boundaryDepth = localDepths(maxIdx);
else
    maxLocalScore = nan;
    boundaryDepth = nan;
end

% 전체 transition score의 robust z-score
scoreMed = median(transitionScore, "omitmissing");
scoreMAD = mad(transitionScore, 1);
boundaryZ = (maxLocalScore - scoreMed) / (scoreMAD + eps);

hasBoundary = isfinite(boundaryZ) && boundaryZ >= changeZThreshold;

fprintf("\n--- Boundary detection ---\n");
fprintf("Expected boundary: %.1f um\n", expectedBoundaryUm);
fprintf("Detected boundary: %.1f um\n", boundaryDepth);
fprintf("Boundary z-score: %.2f\n", boundaryZ);
fprintf("Has reliable boundary? %d\n", hasBoundary);

% 상층 depth masking 생성
keepDepth = true(size(depth_in_um));

if hasBoundary
    switch excludeSide
        case "shallower"
            keepDepth = depth_in_um >= boundaryDepth;
        case "deeper"
            keepDepth = depth_in_um <= boundaryDepth;
        otherwise
            error("excludeSide must be either 'shallower' or 'deeper'.");
    end
else
    warning("Reliable boundary was not detected. No depth exclusion applied.");
end

excludedDepth = ~keepDepth;

fprintf("Excluded depths: %d / %d\n", sum(excludedDepth), numel(excludedDepth));

V_depth_masked = V_depth(keepDepth, :);
depth_masked_in_um = depth_in_um(keepDepth);

%% 11-1. QC figure1: frequency-band depth profile

% figure;
% 
% imagesc(1:nBand, depth_in_um, logMeanBandPSD);
% axis xy;
% xticks(1:nBand);
% xticklabels(bandNames);
% xtickangle(45);
% ylabel("Depth (\mum)");
% title("Depth-wise mean PSD by frequency band");
% colorbar;
% 
% if hasBoundary
%     hold on;
%     yline(boundaryDepth, "w--", "LineWidth", 2);
% end

%% 11-2. transition score
% figure;
% 
% plot(transitionScore, midDepth, "k-", "LineWidth", 1.5);
% set(gca, "YDir", "normal");
% xlabel("Transition score");
% ylabel("Depth between adjacent contacts (\mum)");
% title("Depth-wise transition score");
% 
% if hasBoundary
%     hold on;
%     yline(boundaryDepth, "r--", "LineWidth", 2);
% end

%% 12. depth 방향 smoothing
% boundary가 detection됨 -> masked depth만 사용
% boundary가 detection되지 않음 -> 원본 depth 사용
V_for_csd = V_depth_masked;
depth_for_csd_um = depth_masked_in_um;

if numel(depth_for_csd_um) < 5
    error("Too few depth points remain after masking. Check boundary detection or excludeSide.");
end

V_depth_smooth = smoothdata(V_for_csd, 1, "gaussian", 5);

%% 13. CSD calculation
dz_in_um = median(diff(depth_for_csd_um));

if any(abs(diff(depth_for_csd_um) - dz_in_um) > 1e-6)
    warning("Depth spacing is not perfectly uniform. Median dz is used.");
end

csd = nan(size(V_depth_smooth));

for d = 2:size(V_depth_smooth,1)-1
    csd(d,:) = (V_depth_smooth(d+1,:)-2*V_depth_smooth(d,:)+V_depth_smooth(d-1,:)) / dz_in_um^2;
end

%% 14. CSD smoothing
csd_smooth = smoothdata(csd, 1, "gaussian", 5);
csd_smooth = smoothdata(csd_smooth, 2, "gaussian", 7);

%% 15. heatmap 그리기
t = (0:size(csd_smooth,2)-1) / fs; % 개수 dimension을 sampling rate으로 나누어 다시 sec 차원으로 변환

if opts.makeFigure
figure;
imagesc(t, depth_for_csd_um, csd_smooth);
axis xy;

xlabel("Time (s)");
ylabel("Depth (\mum)");
title("CSD heatmap after superficial layer masking", "Interpreter", "none");
colorbar;

% color scaling: NaN 제외하고 99 percentile 기준
validCsd = csd_smooth(isfinite(csd_smooth));

if ~isempty(validCsd)
    lim = prctile(abs(validCsd), 99);

    if isfinite(lim) && lim > 0
        clim([-lim lim]);
    end
else
    warning("csd_smooth contains no finite values. Color scaling was not applied.");
end

% boundary line 표시: 검출되었고, 현재 plot depth 범위 안에 있을 때만
if hasBoundary && isfinite(boundaryDepth)
    hold on;

    if boundaryDepth >= min(depth_for_csd_um) && boundaryDepth <= max(depth_for_csd_um)
        yline(boundaryDepth, "w--", "LineWidth", 2);
    end
end
end

out.csd_smooth = csd_smooth;
out.fs = fs;
out.t = (0:size(csd_smooth,2)-1) / fs;
out.depth_for_csd_um = depth_for_csd_um;
out.startSec = opts.startSec;
out.durationSec = durationSec;
out.fullDurationSec = fullDurationSec;

end
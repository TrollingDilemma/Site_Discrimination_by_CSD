function result = SCvsRSC_single_from_out(out, weights, makeFigure, verbose)
%SCvsRSC_SINGLE_FROM_OUT Compute one-switch RSC/SC boundary from one CSD segment.
% result.boundaryDepth_oneSwitch is the key output.

if nargin < 2; weights = []; end
if nargin < 3 || isempty(makeFigure); makeFigure = false; end
if nargin < 4 || isempty(verbose); verbose = true; end

csd_smooth = out.csd_smooth;
fs  = out.fs;
t   = out.t;
depth_um = out.depth_for_csd_um;
CSD = csd_smooth;

%% 2. Feature matrix formation in Z-scale
% 2-1. Amplitude feature
nDepth = size(CSD, 1);
nTime = size(CSD, 2); % CSD_Calculating.m의 durationSec or nSamp 참조

csdRMS = rms(CSD,2,"omitmissing");
csdMeanAbs = mean(abs(CSD), 2, "omitmissing");
csdPeaktoPeak = max(CSD, [], 2, "omitmissing") - min(CSD, [], 2, "omitmissing");

% 2-2. Current source/sink balance feature: 부호 주의
sinkComponent = max(-CSD, 0);
sourceComponent = max(CSD, 0);

sinkPower = sum(sinkComponent, 2, "omitmissing");
sourcePower = sum(sourceComponent, 2, "omitmissing");

meanSink = mean(sinkComponent, 2, "omitmissing");
meanSource = mean(sourceComponent,2 , "omitmissing");

sinkSourceRatio = sinkPower ./ (sourcePower + eps);
sinkSourceBalance = (sourcePower - sinkPower) ./ (sourcePower + sinkPower + eps);

% 2-3. Temporal variability feature
csdTemporalSD = std(CSD, 0 ,2, "omitmissing");
csdTemporalMAD = mad(CSD, 1, 2);
dCSDdt = diff(CSD, 1, 2) * fs;
csdTemporalDerivativeRMS = rms(dCSDdt, 2, "omitmissing");

% 2-4. Frequency domain feature
freqBands = struct();
freqBands.delta = [1 4];
freqBands.theta = [4 8];
freqBands.alpha_beta = [8 30];
freqBands.lowGamma = [30 80];
freqBands.highGamma = [80 150];

bandNames = fieldnames(freqBands);
nBand = numel(bandNames);

winSec = 2;
winLength = round(winSec * fs);
winLength = min(winLength, nTime);

noverlap = round(0.5 * winLength);
nfft = max(2^10, 2^nextpow2(winLength));

for d = 1:nDepth
    x = CSD(d,:);
    x = double(x);

    valid = isfinite(x);
    x = x(valid);

    if numel(x) < 10
        continue;
    end

    x = x - mean(x, "omitmissing");

    [Pxx, F] = pwelch(x, winLength, noverlap, nfft, fs);

    for b = 1:nBand
        band = freqBands.(bandNames{b});
        idx = F >= band(1) & F <= band(2);

        if any(idx)
            csdMeanBandPSD(d, b) =mean(Pxx(idx), "omitmissing");
            csdBandPower(d, b) = trapz(F(idx), Pxx(idx));
        end
    end
end

logCsdMeanBandPSD = 10 * log10(csdMeanBandPSD + eps);
logCsdBandPower = 10* log10(csdBandPower + eps);

% 2-5. Adjacent-depth similarity feature
adjCorr = nan(nDepth - 1, 1);

for d = 1:nDepth-1
    x1 = CSD(d,:);
    x2 = CSD(d+1,:);

    valid = isfinite(x1) & isfinite(x2);
    if sum(valid) > 10
        adjCorr(d) = corr(x1(valid).', x2(valid).', "Rows", "complete");
    end
end

corrDrop = 1 - adjCorr;

adjCorrPerDepth = nan(nDepth, 1);
corrDropPerDepth = nan(nDepth, 1);

adjCorrPerDepth(1) = adjCorr(1);
adjCorrPerDepth(end) = adjCorr(end);

corrDropPerDepth(1) = corrDrop(1);
corrDropPerDepth(end) = corrDrop(end);

for d = 2:nDepth-1
    adjCorrPerDepth(d) = mean([adjCorr(d-1), adjCorr(d)], "omitmissing");
    corrDropPerDepth(d) = mean([corrDrop(d-1), corrDrop(d)], "omitmissing");
end

% Defining feature matrix
featureMat = [
    csdRMS, ...
    csdMeanAbs, ...
    csdPeaktoPeak, ...
    sinkPower, ...
    sourcePower, ...
    meanSink, ...
    meanSource, ...
    sinkSourceRatio, ...
    sinkSourceBalance, ...
    csdTemporalSD, ...
    csdTemporalMAD, ...
    csdTemporalDerivativeRMS, ...
    logCsdMeanBandPSD, ...
    logCsdBandPower, ...
    adjCorrPerDepth, ...
    corrDropPerDepth, ...
] ;

featureNames = [
    "csdRMS", ...
    "csdMeanAbs", ...
    "csdPeaktoPeak", ...
    "sinkPower", ...
    "sourcePower", ...
    "meanSink", ...
    "meanSource", ...
    "sinkSourceRatio", ...
    "sinkSourceBalance", ...
    "csdTemporalSD", ...
    "csdTemporalMAD", ...
    "csdTemporalDerivativeRMS", ...
] ;

for b = 1:nBand
    featureNames(end+1) = "meanPSD_" + string(bandNames{b});
end

for b = 1:nBand
    featureNames(end+1) = "bandPower_" + string(bandNames{b});
end

featureNames(end+1) = "adjCorrPerDepth" ;
featureNames(end+1) = "corrDropPerDepth" ;

% Feature matrix의 원소들을 column 방향으로 robust-Z normalization
featureMed = median(featureMat, 1, "omitmissing");
featureMAD = mad(featureMat, 1, 1);

featureMAD(featureMAD == 0 | ~isfinite(featureMAD)) = eps;

featureZ = (featureMat - featureMed) ./ featureMAD;
featureZ(~isfinite(featureZ)) = NaN;

% QC
% fprintf("\n--- CSD feature extraction complete ---\n");
% fprintf("CSD size: %d depths x %d timepoints\n", nDepth, nTime);
% fprintf("Feature matrix size: %d depths x %d features\n", ...
%     size(featureMat,1), size(featureMat,2));
% 
% disp("Feature names:");
% disp(featureNames.');

%% 3. Transition score
% 3-1. Adjacent-depth feature gradient (L2 norm)

featureZ_filled = featureZ; % Nan에 의한 계산 중지 예방 코드

for k = 1:size(featureZ_filled, 2)
    x = featureZ_filled(:, k);
    medVal = median(x, "omitmissing");

    if ~isfinite(medVal)
        medVal = 0;
    end

    x(~isfinite(x))= medVal; % Nan이 있는 feature는 column median으로 대체
    featureZ_filled(:, k) = x;
end

featureDiff = diff(featureZ_filled, 1, 1); % 이웃 depth 별 feature 차분
featureGradient = sqrt(sum(featureDiff.^2, 2)); % L2 norm
midDepth_um= (depth_um(1:end-1) + depth_um(2:end)) / 2; % d_th, (d+1)_th depth의 midpoint

% 3-2. Adjacent-depth CSD correlation drop
% 앞의 corrDrop 사용

% 3-3. Normalizaiton of adjacent-depth feature gradient
fgMed = median(featureGradient,"omitmissing");
fgMAD = mad(featureGradient, 1);

if fgMAD == 0 || ~isfinite(fgMAD)
    fgMAD = eps;
end

fgZ = (featureGradient - fgMed) ./ fgMAD;

% 3-4. Normalization of adjacent-depth CSD correlation drop
cdMed = median(corrDrop, "omitmissing");
cdMAD = mad(corrDrop, 1);

if cdMAD == 0 || ~isfinite(cdMAD)
    cdMAD = eps;
end

corrDropZ = (corrDrop - cdMed) ./ cdMAD;

fgZ(~isfinite(fgZ)) = 0; % Nan이 있는 element는 0으로 대체
corrDropZ(~isfinite(corrDropZ)) = 0; % Nan이 있는 element는 0으로 대체

% 3-5. Summation of two Z's
transitionScore = fgZ + corrDropZ;

%% 4. Multiple boundary candidate detection & segmentation
% 4-1. finding boundary candidates
boundaryZThreshold = 2.5;

tsMed = median(transitionScore, "omitmissing");
tsMAD = mad(transitionScore, 1);

if tsMAD == 0 || ~isfinite(tsMAD)
    tsMAD = eps;
end

tsZ = (transitionScore - tsMed) ./ tsMAD;
tsZ(~isfinite(tsZ)) = 0;

boundaryIdxAll = find(tsZ >= boundaryZThreshold);
boundaryDepthAll = midDepth_um(boundaryIdxAll);
boundaryScoreAll = transitionScore(boundaryIdxAll);
boundaryZAll = tsZ(boundaryIdxAll);

% QC script
% fprintf("\n--- Multiple boundary candidate detection ---\n");
% fprintf("Boundary Z threshold: %.2f\n", boundaryZThreshold);
% fprintf("Number of boundary candidates: %d\n", numel(boundaryIdxAll));
% 
% if ~isempty(boundaryIdxAll)
%     boundaryTable = table( ...
%         boundaryIdxAll(:), ...
%         boundaryDepthAll(:), ...
%         boundaryScoreAll(:), ...
%         boundaryZAll(:), ...
%         'VariableNames', ...
%         {'BoundaryIndex', 'BoundaryDepth_um', 'TransitionScore', 'TransitionScoreZ'} ...
%         );
% 
%     disp(boundaryTable);
% else
%     warning("No boundary candidate exceeded the Z threshold.");
% end

% 4-2. Merging adjacent boundary candidates
if isempty(boundaryIdxAll)
    boundaryZoneStartIdx = [];
    boundaryZoneEndIdx = [];
    boundaryZonePeakIdx = [];
    boundaryZonePeakDepth = [];
    boundaryZonePeakZ = [];
    boundaryZonePeakScore = [];
else
    % boundaryIdxAll에서 연속 구간 찾기
    gap = diff(boundaryIdxAll);

    zoneStartPos = [1; find(gap > 1) + 1];
    zoneEndPos   = [find(gap > 1); numel(boundaryIdxAll)];

    nZone = numel(zoneStartPos);

    boundaryZoneStartIdx = nan(nZone, 1);
    boundaryZoneEndIdx   = nan(nZone, 1);
    boundaryZonePeakIdx  = nan(nZone, 1);
    boundaryZonePeakDepth = nan(nZone, 1);
    boundaryZonePeakZ = nan(nZone, 1);
    boundaryZonePeakScore = nan(nZone, 1);

    for z = 1:nZone
        thisIdx = boundaryIdxAll(zoneStartPos(z):zoneEndPos(z));

        boundaryZoneStartIdx(z) = thisIdx(1);
        boundaryZoneEndIdx(z)   = thisIdx(end);

        % 해당 boundary zone 내부에서 tsZ가 최대인 지점 선택
        [~, localMaxPos] = max(tsZ(thisIdx));
        peakIdx = thisIdx(localMaxPos);

        boundaryZonePeakIdx(z) = peakIdx;
        boundaryZonePeakDepth(z) = midDepth_um(peakIdx);
        boundaryZonePeakZ(z) = tsZ(peakIdx);
        boundaryZonePeakScore(z) = transitionScore(peakIdx);
    end

    boundaryZoneTable = table( ...
        boundaryZoneStartIdx, ...
        boundaryZoneEndIdx, ...
        midDepth_um(boundaryZoneStartIdx), ...
        midDepth_um(boundaryZoneEndIdx), ...
        boundaryZonePeakIdx, ...
        boundaryZonePeakDepth, ...
        boundaryZonePeakScore, ...
        boundaryZonePeakZ, ...
        'VariableNames', ...
        {'StartIndex', 'EndIndex', ...
        'StartDepth_um', 'EndDepth_um', ...
        'PeakIndex', 'PeakDepth_um', ...
        'PeakTransitionScore', 'PeakTransitionScoreZ'} ...
        );
    % fprintf("\n--- Merged boundary zones ---\n");
    % disp(boundaryZoneTable);
end

% 4-3. Initial clustering
boundaryDepthsForSplit = boundaryZonePeakDepth;
boundaryDepthsForSplit = sort(boundaryDepthsForSplit(:));

initialClusterLabel = ones(nDepth, 1);

for b = 1:numel(boundaryDepthsForSplit)
    initialClusterLabel(depth_um > boundaryDepthsForSplit(b)) = b + 1;
end % multiple boundaries에 대해 cluster 1 -> 2 -> 3 ...으로 분류

fprintf("\n--- Initial cluster split ---\n");
fprintf("Number of boundary zones: %d\n", numel(boundaryDepthsForSplit));
fprintf("Number of initial clusters: %d\n", numel(boundaryDepthsForSplit) + 1);

% 4-4. Peri-boundary zone을 ambiguous zone으로 설정
if ~isempty(boundaryZonePeakIdx)
    for b = 1:numel(boundaryZonePeakIdx)
        peakIdx = boundaryZonePeakIdx(b);
        ambiguousIdx = max(1, peakIdx) : min(nDepth, peakIdx + 1);
        initialClusterLabel(ambiguousIdx) = 0; % boundary zone은 cluster 0으로 정의
    end
end

% QC
% figure;
% 
% scatter(initialClusterLabel, depth_um, 40, initialClusterLabel, "filled");
% set(gca, "YDir", "normal");
% xlabel("Initial cluster label");
% ylabel("Depth (\mum)");
% title("Initial cluster split by multiple transition boundaries");

%% 5.Retrosplenial cortex | Superior colliculus discrimination
% 5-1. 필요한 feature column 찾기
findFeature = @(name) find(strcmp(featureNames, name), 1);

idx_csdRMS = findFeature("csdRMS");
idx_csdMeanAbs = findFeature("csdMeanAbs");
idx_sinkSourceBalance = findFeature("sinkSourceBalance");
idx_sinkSourceRatio = findFeature("sinkSourceRatio");
idx_theta = findFeature("bandPower_theta");
idx_lowGamma = findFeature("bandPower_lowGamma");
idx_highGamma = findFeature("bandPower_highGamma");
idx_adjCorrPerDepth = findFeature("adjCorrPerDepth");
idx_corrDropPerDepth = findFeature("corrDropPerDepth");

% 5-2. Segment profile 일람
segmentIDs = unique(initialClusterLabel);
segmentIDs = segmentIDs(segmentIDs > 0); % 0은 ambiguous zone
nSeg = numel(segmentIDs);

segmentID_col = nan(nSeg, 1);
startDepth_um_col = nan(nSeg, 1);
endDepth_um_col = nan(nSeg, 1);
centerDepth_um_col = nan(nSeg, 1);
segmentWidth_um_col = nan(nSeg, 1);
nDepthInSegment_col = nan(nSeg, 1);

% 아래 7가지 항목들은 생물학적 특성에 대한 값들
segAmplitude_col = nan(nSeg, 1);
segSinkSourceAsym_col = nan(nSeg, 1);
segThetaGamma_col = nan(nSeg, 1);
segInternalCorr_col = nan(nSeg, 1);
segInternalGradient_col = nan(nSeg, 1);
segBoundarySharpness_col = nan(nSeg, 1);
segLocalPeak_col = nan(nSeg, 1);

for s = 1:nSeg
    segID = segmentIDs(s);
    idxDepth = find(initialClusterLabel == segID); % 자료형은 cell array

    segmentID_col(s) = segID;
    startDepth_um_col(s) = min(depth_um(idxDepth));
    endDepth_um_col(s) = max(depth_um(idxDepth));
    centerDepth_um_col(s) = median(depth_um(idxDepth));
    segmentWidth_um_col(s) = max(depth_um(idxDepth)) - min(depth_um(idxDepth));
    nDepthInSegment_col(s) = numel(idxDepth);

    % 5-2-1. Amplitude evidence
    rmsVal = NaN;
    meanAbsVal = NaN;

    if ~isempty(idx_csdRMS)
        rmsVal = median(featureZ(idxDepth, idx_csdRMS), "omitmissing");
    end

    if ~isempty(idx_csdMeanAbs)
        meanAbsVal = median(featureZ(idxDepth, idx_csdMeanAbs), "omitmissing");
    end

    segAmplitude_col(s) = median([rmsVal, meanAbsVal], "omitmissing");

    % 5-2-2. Sink/source asymmetry evidence
    ssBal = NaN;
    ssRat = NaN;

    if ~isempty(idx_sinkSourceBalance)
        ssBal = median(featureZ(idxDepth, idx_sinkSourceBalance), "omitmissing");
    end

    if ~isempty(idx_sinkSourceRatio)
        ssRat = median(featureZ(idxDepth, idx_sinkSourceRatio), "omitmissing");
    end

    segSinkSourceAsym_col(s) = median([abs(ssBal), abs(ssRat)], "omitmissing");

    % 5-2-3. Theta/gamma network-like evidence
    thetaVal = NaN;
    lowGammaVal = NaN;
    highGammaVal = NaN;

    if ~isempty(idx_theta)
        thetaVal = median(featureZ(idxDepth, idx_theta), "omitmissing");
    end

    if ~isempty(idx_lowGamma)
        lowGammaVal = median(featureZ(idxDepth, idx_lowGamma), "omitmissing");
    end

    if ~isempty(idx_highGamma)
        highGammaVal = median(featureZ(idxDepth, idx_highGamma), "omitmissing");
    end

    segThetaGamma_col(s) = median([thetaVal, lowGammaVal, highGammaVal], "omitmissing");

    % 5-2-4. Segment-internal adjacent-depth correlation
    pairIdx = idxDepth(idxDepth < nDepth);
    pairIdx = pairIdx(ismember(pairIdx + 1, idxDepth));

    if ~isempty(pairIdx)
        segInternalCorr_col(s) = median(adjCorr(pairIdx), "omitmissing");
    elseif ~isempty(idx_adjCorrPerDepth)
        segInternalCorr_col(s) = median(featureZ(idxDepth, idx_adjCorrPerDepth), "omitmissing");
    else
        segInternalCorr_col(s) = NaN;
    end

    % 5-2-5. Segment-internal feature gradient
    if ~isempty(pairIdx)
        segInternalGradient_col(s) = median(featureGradient(pairIdx), "omitmissing");
    else
        segInternalGradient_col(s) = NaN;
    end

    % 5-2-6. Boundary sharpness
    startIdx = idxDepth(1);
    endIdx   = idxDepth(end);

    boundaryIdx = [];

    if startIdx > 1
        boundaryIdx(end+1) = startIdx - 1;
    end

    if endIdx < nDepth
        boundaryIdx(end+1) = endIdx;
    end

    boundaryIdx = boundaryIdx(boundaryIdx >= 1 & boundaryIdx <= numel(tsZ));

    if ~isempty(boundaryIdx)
        segBoundarySharpness_col(s) = median(tsZ(boundaryIdx), "omitmissing");
    else
        segBoundarySharpness_col(s) = NaN;
    end

    % 5-2-7. Local peak evidence: segment 내부에서 amplitude가 국소적으로 튀는 정도
    if ~isempty(idx_csdRMS)
        ampDepth = featureZ(idxDepth, idx_csdRMS);
        segLocalPeak_col(s) = max(ampDepth, [], "omitmissing") - ...
            median(ampDepth, "omitmissing");
    else
        segLocalPeak_col(s) = NaN;
    end
end

% 5-3. Segment profile matrix formation
segmentProfileMat = [
    segmentWidth_um_col, ...
    segInternalCorr_col, ...
    segInternalGradient_col, ...
    segThetaGamma_col, ...
    segAmplitude_col, ...
    segSinkSourceAsym_col, ...
    segBoundarySharpness_col, ...
    segLocalPeak_col ...
    ];

segmentProfileNames = [
    "segmentWidth_um", ...
    "internalCorr", ...
    "internalGradient", ...
    "thetaGamma", ...
    "amplitude", ...
    "sinkSourceAsymmetry", ...
    "boundarySharpness", ...
    "localPeak" ...
    ];

segMed = median(segmentProfileMat, 1, "omitmissing");
segMAD = mad(segmentProfileMat, 1, 1);

segMAD(segMAD == 0 | ~isfinite(segMAD)) = eps;

spZ = (segmentProfileMat - segMed) ./ segMAD;
spZ(~isfinite(spZ)) = 0; % segment-wise robust-Z normalization

% 5-4. Depth prior
% RSC가 shallower, SC가 deeper라고 가정할 때만 사용
% Depth가 얕을 수록 RSC score에 가중을 주고, 깊을 수록 SC에 가중을 줌

useDepthPrior = true;
depthPriorWeight = 1.5;

if useDepthPrior

    depthNorm = (centerDepth_um_col - min(centerDepth_um_col)) ./ ...
                (max(centerDepth_um_col) - min(centerDepth_um_col) + eps);

    % shallow = RSC-like prior
    superficialPrior = 1 - depthNorm;

    % deep = SC-like prior
    deepPrior = depthNorm;

    superficialPriorZ = (superficialPrior - median(superficialPrior, "omitmissing")) ./ ...
                        (mad(superficialPrior, 1) + eps);

    deepPriorZ = (deepPrior - median(deepPrior, "omitmissing")) ./ ...
                 (mad(deepPrior, 1) + eps);

else
    superficialPriorZ = zeros(nSeg, 1);
    deepPriorZ = zeros(nSeg, 1);
end
% 5-5. Defining RSC score & SC score
    % segmentFeatureZ column:
    % 1 = width
    % 2 = internalCorr
    % 3 = internalGradient
    % 4 = thetaGamma
    % 5 = amplitude
    % 6 = sinkSourceAsymmetry
    % 7 = boundarySharpness
    % 8 = localPeak
    % Weight values는 실제 data와 correlation하며 수정 요망

% weights fields:
% RSC = [width, internalCorr, internalGradient, thetaGamma, boundarySharpness]
% SC  = [amplitude, sinkSourceAsymmetry, localPeak, internalGradient]
% depthPrior = scalar
if nargin < 2 || isempty(weights)
    weights.RSC = [1.00, 1.00, -1.00, 0.70, -0.50];
    weights.SC  = [1.00, 1.00, 0.70, 0.50];
    weights.depthPrior = 1.50;
end

depthPriorWeight = weights.depthPrior;

RSCscore = ...
    + weights.RSC(1) * spZ(:, 1) ...   % broad segment
    + weights.RSC(2) * spZ(:, 2) ...   % high internal continuity
    + weights.RSC(3) * spZ(:, 3) ...   % smooth depth-wise gradient: usually negative
    + weights.RSC(4) * spZ(:, 4) ...   % theta/gamma network-like activity
    + weights.RSC(5) * spZ(:, 7) ...   % not too boundary-sharp: usually negative
    + depthPriorWeight * superficialPriorZ;

SCscore = ...
    + weights.SC(1) * spZ(:, 5) ...   % strong CSD amplitude
    + weights.SC(2) * spZ(:, 6) ...   % sink/source asymmetry
    + weights.SC(3) * spZ(:, 8) ...   % local peak
    + weights.SC(4) * spZ(:, 3) ...   % abrupt internal feature change
    + depthPriorWeight * deepPriorZ;

scoreDiff = RSCscore - SCscore;
smoothScoreDiff = smoothdata(scoreDiff, "gaussian", 3);

splitScore = nan(nSeg-1, 1);
minSegPerRegion = 2;

for k = 1:nSeg-1

    if k < minSegPerRegion || (nSeg-k) < minSegPerRegion
        continue;
    end % 후보 boundary가 너무 위쪽이나 아래쪽에 생겨서 RSC 영역 혹은 SC 영역이 너무 좁아지는 것을 예방

    leftIdx = 1:k;
    rightIdx = k+1:nSeg;

    rscFit = median(smoothScoreDiff(leftIdx), "omitmissing"); 
    scFit = -median(smoothScoreDiff(rightIdx), "omitmissing");

    splitScore(k) = rscFit + scFit;
end

% electrode의 삽입 과정을 고려하여, labeling이 "one-switch" 방식으로 이루어지도록 설계
if all(~isfinite(splitScore))
    kStar = max(1, floor(nSeg/2));
else
    [~, kStar] = max(splitScore);
end

assignedLabel = strings(nSeg, 1);
assignedLabel(1:kStar) = "RSC-like";
assignedLabel(kStar+1:end) = "SC-like";

boundaryDepth_oneSwitch = mean([centerDepth_um_col(kStar),centerDepth_um_col(kStar+1)]);

segmentLabelPerDepth = strings(nDepth, 1);
segmentLabelPerDepth(:) = "unassigned";

% Depth-wise biological mapping

for s = 1:nSeg
    segID = segmentID_col(s);
    idxDepth = initialClusterLabel == segID;
    segmentLabelPerDepth(idxDepth) = assignedLabel(s);
end

% 기존 multiple boundary에서 생긴 initialClusterLabel == 0은 일단 nearest side로 흡수
idxZero = initialClusterLabel == 0;

segmentLabelPerDepth(idxZero & depth_um <= boundaryDepth_oneSwitch) = "RSC-like";
segmentLabelPerDepth(idxZero & depth_um > boundaryDepth_oneSwitch) = "SC-like";

% True ambiguous를 one-switch boundary 주변 좁은 범위만 지정
ambiguousWidthUm = 25;

idxAmbiguous = abs(depth_um - boundaryDepth_oneSwitch) <= ambiguousWidthUm;
segmentLabelPerDepth(idxAmbiguous) = "ambiguous/boundary";

histologicSearchStart_um = boundaryDepth_oneSwitch - ambiguousWidthUm;
histologicSearchEnd_um   = boundaryDepth_oneSwitch + ambiguousWidthUm;

fprintf("\n--- 조직학적 경계 탐색 범위 ---\n");
fprintf("RSC와 SC의 histologic boundary를 다음 depth 범위에서 찾으십시오:\n");
fprintf("    %.1f um ~ %.1f um\n", ...
    histologicSearchStart_um, histologicSearchEnd_um);
fprintf("RSC와 SC 경계의 중심은 %.1f um 깊이에 위치하고, %.1f um의 폭을 갖습니다.\n", ...
    boundaryDepth_oneSwitch, 2 * ambiguousWidthUm);

result.boundaryDepth_oneSwitch = boundaryDepth_oneSwitch;
result.histologicSearchStart_um = histologicSearchStart_um;
result.histologicSearchEnd_um = histologicSearchEnd_um;
result.segmentLabelPerDepth = segmentLabelPerDepth;
result.assignedLabel = assignedLabel;
result.centerDepth_um_col = centerDepth_um_col;
result.splitScore = splitScore;
result.kStar = kStar;
result.RSCscore = RSCscore;
result.SCscore = SCscore;
result.scoreDiff = scoreDiff;
result.depth_um = depth_um;
result.weights = weights;

if makeFigure
    figure;
    labelNumeric = zeros(nDepth, 1);
    labelNumeric(segmentLabelPerDepth == "RSC-like") = 1;
    labelNumeric(segmentLabelPerDepth == "SC-like") = 2;
    labelNumeric(segmentLabelPerDepth == "ambiguous/boundary") = 0;

    scatter(labelNumeric, depth_um, 40, labelNumeric, "filled");
    set(gca, "YDir", "normal");
    xlabel("Biological label: 1=RSC-like, 2=SC-like, 0=ambiguous/boundary");
    ylabel("Depth (\mum)");
    title("Depth-wise biological labeling");
end

end

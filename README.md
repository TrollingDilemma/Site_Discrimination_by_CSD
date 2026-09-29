# CSD_demo

Neuropixels LFP 데이터에서 깊이별 CSD 관련 신호를 계산하고, RSC(retrosplenial cortex)와 SC(superior colliculus)의 경계 후보를 탐색하는 MATLAB 분석 코드입니다.

여러 시간 구간에서 추정 경계가 일관되게 나오도록 점수 가중치를 탐색합니다. 조직학적 정답으로 학습하거나 검증된 분류기는 아니므로, `RSC-like`와 `SC-like` 출력은 탐색적 라벨로 해석해야 합니다.

## 파일 구성

| 파일 | 역할 |
| --- | --- |
| `CSD_Calculating_batchable.m` | 지정 시간 구간 읽기, LFP 전처리, 깊이별 신호 구성, 초기 마스킹 및 CSD 계산 |
| `SCvsRSC_single_from_out.m` | 한 구간의 CSD 특징으로 RSC-like → SC-like 단일 전환 경계 추정 |
| `SCvsRSC_Discriminating_batchOptimize.m` | 무작위 시간 구간별 경계의 일관성을 높이는 가중치 탐색 및 결과 저장 |
| `SGLX_readMeta.m` | SpikeGLX 메타데이터·바이너리 읽기, gain 보정 및 디지털 채널 처리 클래스 |
| `DemoReadSGLXData.m` | 선택한 바이너리의 처음 1초를 읽고 첫 번째 아날로그 채널을 표시하는 예제 |
| `ks25_phy_toBinary.m` | Kilosort 2.5/Phy 처리 데이터를 unwhitening하여 바이너리·메타데이터로 저장하는 별도 도구 |
| `260415_3_g1_tcat.imec0.lf.bin` | 로컬 분석 대상 신호 데이터 |
| `260415_3_g1_tcat.imec0.lf.meta` | 대응하는 기록 설정, gain 및 전극 좌표 정보 |
| `SCvsRSC_batchOptimization_result.mat` | 배치 최적화 결과 파일 |

바이너리와 결과 파일은 코드와 별도로 보관할 수 있으며, GitHub 저장소에 포함되어 있다고 가정하지 않습니다.

## 현재 예제 데이터

로컬 `.meta` 파일에 기록된 정보입니다.

| 항목 | 값 |
| --- | --- |
| 샘플링 속도 | 2,500 Hz |
| 기록 길이 | 1,941.7088초, 약 32분 22초 |
| 저장 채널 | 385개: 신경 신호 384개 + 동기화 1개 |
| 바이너리 크기 | 3,737,789,440 bytes, 약 3.74 GB |
| 데이터 형식 | 채널 × 시간으로 읽는 int16 데이터 |
| 전극 모델 표기 | NP2004 |
| 전극 y 좌표 범위 | 0–2,865 µm |
| 깊이 간격 | 15 µm, 동일 깊이의 두 채널을 평균 |

파일명은 `.lf.bin`이지만 메타데이터의 `snsApLfSy`는 `384,0,1`입니다. 파일명과 AP/LF 채널 구분 표기가 일치하지 않으므로, 다른 데이터에 적용할 때 실제 저장 설정을 확인해야 합니다. 메타데이터에는 CatGT 처리 이력도 포함되어 있습니다.

## 실행 환경

- MATLAB 및 이 폴더에 포함된 `.m` 파일
- 필터링·스펙트럼 분석 함수: `butter`, `filtfilt`, `pwelch`, `designNotchPeakIIR` 등
- 통계 관련 함수: `mad`, `corr`, `prctile` 등

Signal Processing Toolbox와 Statistics and Machine Learning Toolbox의 관련 함수 사용 가능 여부를 확인하십시오. `designNotchPeakIIR` 및 `"omitmissing"` 등의 문법을 지원하는 환경이 필요합니다. 최소 MATLAB 릴리스와 툴박스 조합은 아직 실행 검증하지 않았습니다.

설치 상태와 주요 함수 경로는 다음과 같이 확인할 수 있습니다.

```matlab
ver
which butter -all
which designNotchPeakIIR -all
which pwelch -all
which mad -all
which corr -all
```

`ks25_phy_toBinary`는 주 분석 흐름에서 사용하지 않습니다. 별도 사용 시 Phy 출력 폴더, `npy-matlab`의 `readNPY`/`writeNPY`, 올바른 unwhitening을 위한 `rez2.mat`이 필요합니다. C_Waves 실행 옵션은 기본적으로 꺼져 있습니다.

## 빠른 시작

### 1. 경로와 데이터 준비

MATLAB의 현재 폴더를 프로젝트 폴더로 설정하고, `.bin`과 같은 이름의 `.meta` 파일을 함께 준비합니다.

```matlab
cd('C:\Users\DohyunShin\Documents\MATLAB\CSD_demo');
```

**현재 `CSD_Calculating_batchable.m`은 함수 인자로 받은 `binPath`를 내부 절대경로로 덮어씁니다.** 다른 위치나 데이터에서 실행하려면 해당 함수의 `%% 1. 경로와 파일명` 아래에 있는 `binPath = ...` 줄을 수정하거나 제거해야 합니다. `metaPath` 인자는 현재 사용되지 않으며, `.bin`과 같은 폴더에서 같은 기본 이름의 `.meta`를 자동으로 읽습니다.

### 2. 짧은 신호 읽기 확인

```matlab
DemoReadSGLXData
```

파일 선택창에서 바이너리를 선택하면 처음 1초의 첫 번째 채널을 gain 보정하여 µV 단위로 표시합니다.

### 3. 30초 구간의 CSD 및 경계 계산

```matlab
binPath = "260415_3_g1_tcat.imec0.lf.bin";
metaPath = "260415_3_g1_tcat.imec0.lf.meta";

opts = struct();
opts.startSec = 0;
opts.durationSec = 30;
opts.makeFigure = true;

out = CSD_Calculating_batchable(binPath, metaPath, opts);
result = SCvsRSC_single_from_out(out, [], true, true);

disp(result.boundaryDepth_oneSwitch);
disp([result.histologicSearchStart_um, result.histologicSearchEnd_um]);
```

빈 가중치 `[]`는 기본 가중치를 사용합니다. `out.t`는 선택 구간 시작을 0초로 하는 상대 시간이며, 파일 기준 시간은 `out.startSec + out.t`입니다.

`opts.durationSec`를 생략하면 전체 기록 길이가 선택됩니다. 전체 데이터를 한 번에 읽는 대신 짧은 구간부터 실행하십시오. 현재 파일은 double 배열 하나만 약 15 GB이며, 전처리 중 여러 배열이 추가로 필요합니다. `startSec > 0`일 때는 기록 끝을 넘지 않도록 `durationSec`도 지정해야 합니다.

### 4. 배치 최적화

```matlab
SCvsRSC_Discriminating_batchOptimize
```

이 스크립트는 시작 시 `clear; clc; close all;`을 실행합니다. 작업 공간 변수를 지우고 열린 figure를 닫으므로 필요한 작업을 먼저 저장하십시오.

기본 설정은 다음과 같습니다.

| 설정 | 기본값 |
| --- | --- |
| 난수 시드 | `rng(1)` |
| 구간 길이 | 30초 |
| 배치 수 × 배치당 구간 수 | 10 × 10 = 100개 |
| 가중치 후보 수 | 기본값 1개 + 무작위 후보 50개 |
| 경계 집계 단위 | 20 µm |

시작 시점은 실제 전체 기록 길이에서 무작위로 선택하며, 구간 간 중복을 허용합니다. 배치는 그룹 구분이며, 병렬 실행이나 순차 학습을 뜻하지 않습니다. 100개 CSD 결과를 메모리에 보관한 뒤 모든 후보를 평가하므로 시간과 메모리 사용량이 클 수 있습니다.

완료하면 현재 폴더의 `SCvsRSC_batchOptimization_result.mat`에 결과를 저장합니다. **같은 이름의 기존 결과는 덮어씁니다.** 보존하려면 실행 전에 기존 결과를 다른 이름으로 복사하거나 스크립트의 저장 파일명을 바꾸십시오.

### 5. 저장 결과 확인 및 가중치 재사용

```matlab
whos('-file', 'SCvsRSC_batchOptimization_result.mat')
saved = load('SCvsRSC_batchOptimization_result.mat');

disp(saved.bestWeights);
disp(saved.finalBoundaryDepth_oneSwitch);
disp(saved.batchSummary);

% 앞서 계산한 out에 최적 가중치 적용
result = SCvsRSC_single_from_out(out, saved.bestWeights, true, true);
```

## 분석 방법

### LFP 전처리와 초기 깊이 마스킹

`CSD_Calculating_batchable`의 처리 순서는 다음과 같습니다.

1. 시간 구간 읽기 및 gain 보정: 신경 신호를 µV로 변환합니다.
2. 선형 추세 제거 및 1–300 Hz bandpass를 적용합니다.
3. 60, 120, 180, 240 Hz notch를 적용합니다.
4. 채널 표준편차가 중앙값의 5배 초과 또는 0.2배 미만이면 불량 채널로 판정합니다.
5. 정상 채널이 2개 이상이면 불량 채널을 채널 인덱스 기준으로 선형 보간합니다.
6. 매 시점의 채널 중앙값을 빼서 공통 reference를 제거합니다.
7. shank 0의 사용 전극을 선택하고 같은 y 좌표의 채널을 평균합니다.
8. RMS, 대역별 PSD, 이웃 깊이 간 상관관계로 변화 점수를 계산합니다.
9. 1,500 ± 400 µm에서 최대 변화점을 찾고, robust 점수가 2.5 이상이면 작은 y 좌표 쪽을 제외합니다. 기준을 넘지 못하면 제외하지 않습니다.

여기서 사용하는 주파수 대역은 delta 1–4 Hz, theta 4–8 Hz, alpha/beta 8–30 Hz, low gamma 30–80 Hz, high gamma 80–150 Hz입니다.

### CSD 계산

남은 깊이별 전압을 Gaussian 방식으로 평활화한 후 다음 값을 계산합니다.

```text
CSD(d, t) = [V(d+1, t) - 2V(d, t) + V(d-1, t)] / dz²
```

이후 깊이·시간 방향으로 추가 평활화합니다. 현재 계산은 전압의 **양의 공간 2차 미분**이며, 단위는 µV/µm²입니다. 전도도와 음의 부호를 곱하는 물리적 CSD 정의가 적용되어 있지 않으므로 절대 전류 밀도로 해석하지 마십시오. 후속 코드의 음수=sink, 양수=source 명칭과 사용할 물리적 부호 정의의 일치 여부도 확인해야 합니다.

### RSC–SC 경계 추정

`SCvsRSC_single_from_out`은 CSD의 진폭, 부호별 성분, 시간 변동, 주파수 및 공간적 연속성 특징으로 깊이 구간을 나눕니다. 각 구간에 가중 점수를 부여하고 RSC-like → SC-like로 한 번만 바뀌는 분할을 선택합니다.

| 점수 | 가중치에 사용되는 특징 순서 | 기본 가중치 |
| --- | --- | --- |
| RSC | 구간 폭, 내부 상관, 내부 변화량, theta/gamma 특징, 경계 급격함 | `[1, 1, -1, 0.7, -0.5]` |
| SC | 진폭, sink/source 비대칭, 국소 피크, 내부 변화량 | `[1, 1, 0.7, 0.5]` |
| 깊이 prior | 작은 좌표 쪽 RSC, 큰 좌표 쪽 SC 선호 | `1.5` |

경계는 선택된 두 인접 구간의 **중심 좌표 사이 중간값**입니다. 경계 ±25 µm는 `ambiguous/boundary`로 표시하며 조직학적 탐색 범위로 출력합니다. 이는 고정된 설정값으로, 통계적 신뢰구간이 아닙니다.

### 가중치 최적화

각 후보의 경계값을 20 µm 단위로 반올림하여 집계하고 다음 목적함수를 최소화합니다.

```text
objective = 1000 × (nUnique - 1)
            + 2 × spreadUm
            + 5 × madUm
            - 500 × convergenceRate
```

- `nUnique`: 반올림한 경계의 서로 다른 값 개수
- `spreadUm`: 유효한 원래 경계값의 최댓값−최솟값
- `madUm`: 원래 경계값의 중앙절대편차
- `convergenceRate`: 유효한 경계 중 가장 빈번한 반올림 경계가 차지하는 비율

최종 경계는 최적 후보의 반올림 경계값 중 최빈값입니다. 이 목적함수는 **시간 구간 간 일관성**을 평가하며 해부학적 정확도를 직접 평가하지 않습니다.

## 주요 출력

| 출력 | 내용 |
| --- | --- |
| `out.csd_smooth` | 마스킹·평활화 후 깊이 × 시간 CSD 행렬 |
| `out.depth_for_csd_um` | CSD 행렬 각 행에 대응하는 전극 y 좌표 |
| `out.fs`, `out.t` | 샘플링 속도와 구간 내 상대 시간 |
| `out.startSec`, `out.durationSec`, `out.fullDurationSec` | 선택 구간 및 전체 기록 길이 |
| `result.boundaryDepth_oneSwitch` | 단일 전환 경계 후보 |
| `result.segmentLabelPerDepth` | 깊이별 RSC-like, SC-like, ambiguous/boundary 라벨 |
| `result.histologicSearchStart_um`, `result.histologicSearchEnd_um` | 경계 ±25 µm 범위 |
| `result.RSCscore`, `result.SCscore`, `result.splitScore` | 구간별 점수와 분할 후보 점수 |

배치 저장 코드에는 다음 변수가 포함됩니다: `bestWeights`, `bestBoundaryVec`, `finalBoundaryDepth_oneSwitch`, `resultTable`, `batchSummary`, `optTable`, `startSecMat`, `sampleSec`, `boundaryTolUm`. 전체 `csdOut` 배열은 저장하지 않습니다. 기존 `.mat` 파일의 실제 수치는 이 문서 작성 과정에서 검증하지 않았습니다.

## 현재 구현의 한계

- 전극 y 좌표를 조직 표면 기준 깊이로 변환하지 않습니다. 작은 좌표가 실제로 superficial 쪽인지 삽입 방향과 전극 배치로 확인해야 합니다.
- CSD 전 단계의 1,500 µm 주변 마스킹과 CSD 후 단계의 RSC–SC 경계 추정은 서로 다른 처리입니다. 초기 마스킹으로 제외된 영역은 후속 분류에 포함되지 않습니다.
- 불량 채널 보간은 실제 전극 거리 대신 채널 인덱스를 사용합니다. 다른 채널 배치나 여러 shank 데이터에 그대로 적용하기에 적합하지 않을 수 있습니다.
- 깊이 간격이 불균일해도 중앙 간격으로 미분하며 경고만 출력합니다.
- `ChannelCountsIM`의 반환 순서는 `[AP, LF, SY]`이지만 CSD 함수는 `[LF, AP, SY]`로 받습니다. AP/LF 출력 표시는 뒤바뀔 수 있습니다. 현재 신경 채널 수 계산은 `nChan - SY`를 사용합니다.
- 분할 후보가 부족하면 중앙 분할로 대체합니다. 구간이 하나뿐인 경우에는 `kStar+1` 참조에서 오류가 발생할 수 있습니다.
- `verbose` 인자는 현재 출력문을 제어하지 않아 `false`로 전달해도 일부 메시지가 출력됩니다.
- 이 README는 코드와 메타데이터의 정적 확인을 바탕으로 작성했습니다. MATLAB 실행 성공, 기존 결과 수치 및 조직학적 타당성을 보증하는 테스트 결과는 포함하지 않습니다.

## GitHub에서 관리하기

`.m`과 `README.md`는 텍스트이므로 일반 Git으로 변경 이력을 관리할 수 있습니다. `.meta`는 내부 경로와 장비 정보를 확인한 후 포함할 수 있습니다. 대용량 `.bin`은 별도 저장소에 보관하고, 데이터 위치·접근 방법을 프로젝트 운영 방식에 맞게 문서화하십시오.

아래는 선택적으로 사용할 `.gitignore` 예시입니다. 이 README 작성만으로 `.gitignore`가 생성되지는 않습니다.

```gitignore
# 대용량 신호 데이터
*.bin
*.dat

# MATLAB 편집기 임시 파일
*.asv

# 재생성 가능한 최적화 결과: 필요에 따라 제외 해제
SCvsRSC_batchOptimization_result.mat
```

결과 `.mat`을 커밋할 경우 입력 데이터 식별 정보와 분석 설정도 함께 기록하는 것이 좋습니다. 외부에서 유래한 유틸리티를 재배포할 때는 원본의 라이선스와 고지 조건을 확인하십시오.

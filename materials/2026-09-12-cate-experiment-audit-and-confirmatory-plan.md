# DataFusionPPI CATE 실험 감사와 확인 실험 계획

> 2026-09-13 supersession: 이 문서는 과거 감사와 설계의 역사 기록으로 보존한다.
> 현재 CATE 구현·실험 범위는 [원고 충실 최소 실험 PRD](./2026-09-12-cate-confirmatory-execution-prd.md)의 활성 본문이 대체한다.
> 아래 M0–M5, 네 사분면, STAR 재구성, 확장 실험은 현재 실행의 선행조건이 아니다.
> 통계 식의 최종 authority는 현재 canonical manuscript다. 아래 handoff 우선순위와 실행 지시는 비활성이다.

작성일: 2026-09-12
상태: 실행 전 계획
통계 specification authority: `2026-09-09-implementation-handoff-JA.md` version 3

## 0. Output Contract

### 0.1 Deliverable

이 문서는 현재 CATE 구현과 결과를 추적 가능하게 감사하고, FusionPPI를 가장
강하면서도 공정하게 검증할 simulation과 real-world 실행계획을 정의한다.

### 0.2 Consumer

용한은 이 문서로 다음 연구·구현 순서를 결정한다. 후속 구현자는 각 milestone의
입력, 산출물, acceptance test, 중단 조건을 실행 계약으로 사용한다.

### 0.3 Form

이 문서는 다음 두 역할만 맡는다.

1. 2026-09-12 현재 CATE evidence의 audit와 status snapshot
2. handoff version 3을 따르는 confirmatory execution roadmap

수식과 통계적 정의가 충돌하면 handoff version 3과 manuscript가 우선한다. 이 문서는
새로운 theorem이나 경쟁하는 formula specification이 아니다.

### 0.4 Acceptance criteria

이 문서가 완료되려면 다음을 모두 포함해야 한다.

- 현재 CATE code, completed run, CSV, plot의 정확한 inventory
- 주장 차원별 `PASS`, `PARTIAL`, `FAIL`과 그 근거
- true CATE risk, reporting score, oracle, regret, coverage의 구분
- $\lambda$와 $\omega$의 기회를 독립적으로 식별하는 DGP
- positive, neutral, adversarial, failure regime
- micro fixture부터 confirmatory simulation까지의 단계와 반복 예산
- STAR, NSW/CPS·PSID, WHI 후보의 estimand와 validation 한계
- testable implementation PRD, 표·그림 PRD, provenance 계약

### 0.5 Evidence

사실 주장은 가능한 한 현재 canonical file과 line에 연결한다. CSV 행 수와 요약값은
해당 raw CSV를 read-only로 집계한 값이다.

### 0.6 Boundaries

이 문서 작성은 실험 실행 승인이 아니다. 이번 작업은 이 파일 하나만 새로 만든다.
code, test, CSV, PNG, manuscript, 기존 handoff, 기존 roadmap을 바꾸지 않는다.
dependency 설치, mirror sync, stage, commit, push도 하지 않는다. 향후 routine 구현은
승인된 milestone 범위에서 수행하되, 통계 protocol 변경과 commit·push는 별도 지시를
받는다.

## 1. Executive verdict

### 1.1 가장 강한 현재 결론

현재 결과는 다음 좁은 주장을 지지한다.

> 공정하게 ridge를 선택한 basis-matched RCT-only learner와 비교하면, 현재의
> common-marginal synthetic cells에서 FusionPPI의 selected rule은 많은 반복에서
> 더 낮은 integrated squared CATE risk를 보였다.

현재 monitoring report는 selected/RCT median risk ratio가 spline에서 trial size가
$100,200,400$일 때 각각 $0.231,0.408,0.625$이고, neural에서는
$0.645,0.754,0.834$라고 보고한다. 전체 반복의 88.2%에서 selected rule의 위험이
RCT-only보다 낮았다
([monitoring report:209](./2026-09-12-implementation-monitoring-report.md#L209)).

### 1.2 현재 지지되지 않는 결론

현재 결과는 다음 넓은 주장을 지지하지 않는다.

> FusionPPI가 일반적인 observational confounding에서 단순 pooled learner보다
> 우월하다.

현재 synthetic DGP에서는 confounding bias가 $X$에 따라 거의 변하지 않는다. source
indicator가 이 상수형 편향을 흡수하므로 domain-indicator pooling이 FusionPPI보다
훨씬 낮은 위험을 보인다
([monitoring report:240](./2026-09-12-implementation-monitoring-report.md#L240),
[monitoring report:257](./2026-09-12-implementation-monitoring-report.md#L257)).
이것은 단순히 반복 수가 부족한 문제가 아니다. 현재 DGP가 두 방법을 판별하지 못하는
설계 문제다.

### 1.3 공식 상태

공식 v3 gate 상태는 `G0`, 결정은 `INCONCLUSIVE`다. 과학적 performance 반복은 8개
cell에서 한 번씩인 총 8회다. CATE-1 neural과 CATE-2는 assertion-only smoke다.
`--stage g1`은 `G1_NOT_IMPLEMENTED`로 fail closed한다
([roadmap:65](./implementation-experiment-roadmap.md#L65),
[roadmap:97](./implementation-experiment-roadmap.md#L97)).

따라서 이후 CATE-1, CATE-2, CATE-3 CSV는 중요한 exploratory evidence이지만 공식
gate를 통과한 confirmatory evidence라고 부르지 않는다.

### 1.4 최종 판정표

| 차원 | 판정 | 의미 |
| --- | --- | --- |
| CATE target과 fusion loss | `PASS` | manuscript의 target-preservation 식은 명확하다. |
| 공정한 RCT-only ridge 비교 | `PASS`, exploratory | neural RCT-only도 같은 validation score로 ridge를 고른다. |
| current synthetic에서 RCT-only 대비 risk 감소 | `PASS`, exploratory | 다수 cell에서 일관된 방향이다. |
| domain-indicator pooling 대비 우월성 | `FAIL`, DGP 판별력 | pooled comparator가 현재 synthetic에서 더 좋다. |
| Theorem 6 bounded check | `PASS`, 실질 증거 약함 | 유효하지만 radius가 매우 느슨하다. |
| Theorem 7와 Corollary 7.1 empirical check | `FAIL`, protocol deviation | 현재 spline은 $\rho=10^{-2}$와 clipping을 사용한다. |
| DRF sieve oracle calculus | `PARTIAL` | 식은 있으나 adaptive MCSE oracle budget이 없다. |
| CATE four-role independence | `FAIL`, protocol deviation | current engine은 3-way split과 shared evaluation을 쓴다. |
| STAR exact source law | `FAIL`, protocol deviation | coarse pattern과 in-sample residual을 사용한다. |
| STAR reporting score | `PASS`, exploratory stress metric | true CATE risk 또는 causal truth는 아니다. |
| CATE-3 ratio experiment | `PASS`, exploratory extension | paired five-route 계산은 있으나 theorem extension은 아니다. |
| manuscript empirical section | `FAIL`, stale | 여전히 CATE experiment가 planned라고 적는다. |

## 2. Scope와 용어

### 2.1 CATE estimand

$X$는 treatment 전에 측정한 공변량이고, $Y(1),Y(0)$은 potential outcomes다. 목표
CATE는 RCT target population에서

$$
\tau(x)=\mathbb E_R\{Y(1)-Y(0)\mid X=x\}
$$

이다. candidate function을 $\zeta(x)$라 하면 true CATE risk는

$$
\mathcal R(\zeta)=\mathbb E_R[\{\zeta(X)-\tau(X)\}^2]
$$

이다
([manuscript 4.tex:11](../manuscript/main/4.tex#L11)). Synthetic experiment에서는
$\tau(X)$를 알기 때문에 independent target-law truth sample로 이 값을 계산할 수 있다.

### 2.2 두 channel

$\lambda$는 RCT outcome regression과 OBS outcome regression을 섞어 RCT-valid
pseudo-outcome 또는 R-loss residual의 조건부 noise를 줄이는 channel이다. $\omega$는
OBS에서 학습한 treatment-effect prediction

$$
\widehat g(x)=\widehat\mu_{O,1}(x)-\widehat\mu_{O,0}(x)
$$

과 큰 OBS covariate sample을 이용해 $X$에 따라 예측 가능한 변동을 줄이는 channel이다
([manuscript 4.tex:21](../manuscript/main/4.tex#L21),
[manuscript 4.tex:184](../manuscript/main/4.tex#L184)).

### 2.3 transport ratio

$P_R^X$와 $P_O^X$는 RCT와 OBS의 covariate law다. 공통 marginal이면
$P_R^X=P_O^X$이고 $r_0(x)=1$이다. 다르면 exact transport ratio는

$$
r_0(x)=\frac{dP_R^X}{dP_O^X}(x)
$$

이다. exact $r_0$는 OBS average를 RCT target으로 옮긴다. fitted ratio는 별도 오차와
drift를 만든다
([manuscript 4.tex:68](../manuscript/main/4.tex#L68),
[manuscript 4.tex:361](../manuscript/main/4.tex#L361)).

### 2.4 서로 다른 평가량

다음 수치는 서로 대체할 수 없다.

1. **True CATE risk:** $\mathcal R(\widehat\zeta)$. Synthetic truth가 있을 때만 직접
   계산한다.
2. **Held-out reporting-score difference:** independent reporting sample에서
   selected와 RCT-only의 validation score 차이를 계산한다. 음수면 해당 score에서
   selected가 더 좋다. STAR에서 이것은 true CATE risk가 아니다.
3. **Grid oracle:** 실제 $\tau$를 사용해 finite candidate grid에서 risk가 가장 작은
   candidate를 고른다. 실행 가능한 estimator가 아니라 ceiling이다.
4. **Sieve-moment oracle:** $A_p,B_p,C_p,D_p$의 population moment로 first-order optimal
   coefficient를 계산한다. finite-grid oracle과 다른 대상이다.
5. **Feasible plug-in oracle:** population moment를 independent Monte Carlo 또는 sample
   moment로 근사한 diagnostic이다. true population oracle과 같지 않다.
6. **Selection regret:**
   $\mathcal R(\widehat\zeta_{selected})-\mathcal R(\widehat\zeta_{grid\ oracle})$다.
7. **Pointwise coverage:** 각 $x$에서 interval이 $\tau(x)$를 포함하는 빈도다. integrated
   risk가 작다고 coverage가 맞는 것은 아니다.

## 3. 현재 implementation inventory

### 3.1 공통 engine

중심 함수는
[`fusion_cate.cate_replication`](../code/fusion_cate.py#L26)이다. 현재 다음을 한다.

- RCT와 OBS nuisance fit, OBS propensity fit, ratio fit
  ([fusion_cate.py:38](../code/fusion_cate.py#L38))
- spline 또는 frozen neural basis와 RCT-nuisance standardization
  ([fusion_cate.py:43](../code/fusion_cate.py#L43))
- DRF와 RF에서 finite grid candidate fit
  ([fusion_cate.py:69](../code/fusion_cate.py#L69))
- validation score selection과 true-risk grid oracle 계산
  ([fusion_cate.py:89](../code/fusion_cate.py#L89))
- RCT-only, $\lambda$-only, $\omega$-only, joint restriction 비교
  ([fusion_cate.py:102](../code/fusion_cate.py#L102))
- DRF spline oracle calculus
  ([fusion_cate.py:143](../code/fusion_cate.py#L143),
  [fusion_cate.py:168](../code/fusion_cate.py#L168))
- experimental grounding과 domain-indicator pooling
  ([fusion_cate.py:155](../code/fusion_cate.py#L155))

현재 engine은 RCT와 OBS 전체 sample을 각각 nuisance, tune, eval의 세 역할로 나눈다
([fusion_cate.py:34](../code/fusion_cate.py#L34)). v3는 nuisance, tuning, selection,
reporting의 네 independent roles를 요구한다
([handoff:817](./2026-09-09-implementation-handoff-JA.md#L817)).

### 3.2 CATE-1 completed runs

driver는 [`exp_cate1.py`](../code/exp_cate1.py#L20)다.

| artifact | rows | 독립 반복 단위 |
| --- | ---: | --- |
| `cate1_spline_replications.csv` | 24,000 | 5 families × 6 configurations × 100 replication IDs |
| `cate1_neural_replications.csv` | 12,000 | 5 families × 6 configurations × 100 replication IDs |
| `cate1_spline_summary.csv` | 2,172 | spline long-row summaries |
| `cate1_neural_summary.csv` | 996 | neural long-row summaries |
| `cate1_spline_decomposition.csv` | 360 | family·cell·sieve·learner·rule 집계 |
| `cate1_neural_decomposition.csv` | 180 | family·cell·learner·rule 집계 |

24,000과 12,000은 독립 반복 수가 아니다. 같은 raw replication에서 sieve, learner,
baseline 행이 여러 개 생긴 long-format row 수다. spline과 neural run 각각의
family-configuration-replication 단위는 3,000개다.

v3는 SCM-B, Gaussian SCM 2, Gaussian SCM 3, STAR real-$X$의 4 families와 family당
6 configurations, 총 24 cells를 규정한다
([handoff:794](./2026-09-09-implementation-handoff-JA.md#L794)). Current driver는
`scm1`, `scm2`, `scm3`, `star`, `bounded`의 5 families를 실행한다. 이것은 useful
exploratory extension이지만 official 24-cell run과 같지 않다.

### 3.3 CATE-2 completed runs

driver는 [`exp_cate2.py`](../code/exp_cate2.py#L33)다.

| artifact | rows | 내용 |
| --- | ---: | --- |
| `cate2_star_real_replications.csv` | 6,400 | 2 outcomes × 2 $n_R$ × 2 $\alpha$ × 2 sieves × 100 reps × 4 result rows |
| `cate2_star_real_summary.csv` | 576 | CATE long-row summaries |
| `cate2_star_real_ate_replications.csv` | 8,000 | 같은 STAR source simulation의 ATE rows |
| `cate2_star_real_ate_summary.csv` | 400 | ATE long-row summaries |
| `cate2_star_real_prediction_variability.csv` | 64 | weighted prediction variability summaries |

현재 source law는 27개 coarse pattern의 같은 $X$ law를 양 source에서 사용한다
([fusion_data.py:130](../code/fusion_data.py#L130)). 그러나 eligibility는 arm당 최소 5명
기준이고
([fusion_data.py:155](../code/fusion_data.py#L155),
[fusion_data.py:179](../code/fusion_data.py#L179)), residual은 whole-cohort in-sample
ridge residual이다
([fusion_data.py:189](../code/fusion_data.py#L189)). 이는 exact tuple,
$e_x\in[0.15,0.85]$, OOF rank, frozen fingerprint를 요구하는 v3와 다르다
([handoff:856](./2026-09-09-implementation-handoff-JA.md#L856)).

### 3.4 CATE-3 completed runs

driver는 [`exp_cate3.py`](../code/exp_cate3.py#L1)다.

| artifact | rows | 독립 raw-DGP 반복 |
| --- | ---: | ---: |
| `cate3_shift_spline_replications.csv` | 24,000 | 1,200 |
| `cate3_shift_neural_replications.csv` | 6,000 | 300 |
| `cate3_shift_spline_summary.csv` | 1,320 | spline long-row summaries |
| `cate3_shift_neural_summary.csv` | 330 | neural long-row summaries |
| `cate3_shift_spline_decomposition.csv` | 360 | 집계 artifact |
| `cate3_shift_neural_decomposition.csv` | 90 | 집계 artifact |

spline은 3 families × 2 shifts × 2 $n_R$ × 100 replications다. neural은 3 families ×
strong shift × $n_R=200$ × 100 replications다. 각 raw replication에서 five ratio routes를
paired evaluation하므로 route 수를 독립 반복 수에 곱하지 않는다. 동일 data pairing은
[exp_cate3.py:82](../code/exp_cate3.py#L82)에 명시돼 있다.

현재 ratio route 구현에는 oracle, unit, classifier, BAL-X, BAL-X+$g$가 모두 있다
([exp_ate2.py:42](../code/exp_ate2.py#L42)). BAL-X+$g$는 OBS fit
$\widehat g$를 balance feature에 실제로 추가한다
([exp_ate2.py:30](../code/exp_ate2.py#L30)). 그러나 기존 통합 shift figure의 route 목록은
oracle, classifier, BAL-X, unit만 사용해 BAL-X+$g$를 누락한다
([make_figures.py:70](../code/make_figures.py#L70)). Confirmatory 표와 그림에는 five routes를
모두 포함해야 한다.

### 3.5 현재 plots와 reports

- `figure2_cate_risk.png`: common-marginal selected/RCT risk
- `figure4_cate_variance_bias.png`: anchor prediction variance와 squared bias
- `2026-09-12-cate3-correct-interpretation.png`: CATE-3 ratio routes와 mechanism
- `2026-09-12-implementation-monitoring-report.md`: 현재 통합 수치 해석
- `2026-09-12-traceability-map.md`: claim, code, test, artifact 연결

manuscript empirical section은 여전히 “CATE experiments are planned”라고 적는다
([manuscript 5.tex:28](../manuscript/main/5.tex#L28)). 현재 results를 manuscript evidence로
사용하려면 confirmatory protocol을 먼저 완료하고 이 불일치를 별도 승인으로 고쳐야 한다.

## 4. 차원별 감사

### 4.1 Target preservation

**판정: `PASS`, theorem statement.**

DRF와 RF loss는 exact transport에서 candidate-dependent target을
$\mathcal R(\zeta)$로 보존한다. $\omega$ correction의 기대값은 0이고 $\lambda$는 noise
floor를 바꾸지만 target CATE를 바꾸지 않는다
([manuscript 4.tex:68](../manuscript/main/4.tex#L68)).

이 계획에서 공식 **current CATE theorem-check 범위는 common marginal,
$r_0=1$**이다. Manuscript의 fixed exact-ratio target identity는 별도로 유효하지만,
current CATE-3는 fitted ratio의 finite-sample behavior를 보여주는 exploratory extension이다.
Estimated-ratio learning equation의 drift를 포함하는 완전한 CATE risk theorem과
confirmatory protocol이 추가되기 전에는 CATE-3를 current common-marginal theorem의 직접
검증이라고 부르지 않는다.

### 4.2 Honest selection and reporting

**판정: `FAIL`, v3 protocol deviation.**

manuscript algorithm은 nuisance, tuning, evaluation을 분리한다
([manuscript 4.tex:158](../manuscript/main/4.tex#L158)). v3는 이를 더 엄격히 하여
selection과 reporting을 독립 역할로 분리한다. Current CATE-1 engine은 하나의 `eval`로
candidate selection과 reported score를 함께 만든다
([fusion_cate.py:56](../code/fusion_cate.py#L56),
[fusion_cate.py:91](../code/fusion_cate.py#L91)). Synthetic true risk는 independent test
truth를 사용하므로 selected-vs-RCT descriptive risk 자체는 계산 가능하지만, selection
score calibration과 unbiased reporting-score claim은 v3 역할 계약을 만족하지 않는다.

### 4.3 Spline theorem check

**Theorem 6 판정: `PASS`, 매우 느슨함.**

SCM-B clipped bounded path만 고정 $B_0=9$, $\bar r_0=1$을 사용한다. Gaussian과 STAR는
radius N/A여야 한다
([handoff:785](./2026-09-09-implementation-handoff-JA.md#L785)). Current report의 bounded
radius는 observed regret의 약 130배여서 유효하지만 substantive performance evidence는
거의 제공하지 않는다
([monitoring report:429](./2026-09-12-implementation-monitoring-report.md#L429)).

**Theorem 7와 Corollary 7.1 판정: `FAIL`, current result label.**

v3 eligible path는 spline, $\rho=0$, unclipped, $p<n_R^{tune}$, condition number
$\le10^8$이다. Current spline ridge grid는 $(10^{-2},)$이고
([fusion_cate.py:43](../code/fusion_cate.py#L43)), candidate와 risk prediction을 clip한다
([fusion_cate.py:73](../code/fusion_cate.py#L73),
[fusion_cate.py:89](../code/fusion_cate.py#L89)). Current rows는 empirical fixed-ridge
candidate 결과로 보존하고 theorem-check table에서 제외해야 한다.

### 4.4 Neural path

**판정: `PASS`, empirical candidate.**

neural representation은 finite grid empirical extension이다. linear-sieve Theorem 7이나
Corollary 7.1의 check가 아니다
([handoff:787](./2026-09-09-implementation-handoff-JA.md#L787),
[manuscript 4.tex:236](../manuscript/main/4.tex#L236)). RCT-only neural candidate가 같은
score로 ridge를 고르는 점은 fair comparison이다
([fusion_cate.py:102](../code/fusion_cate.py#L102)).

### 4.5 Oracle correctness

**판정: `PARTIAL`.**

Current DRF spline code는 $A_p,B_p,C_p,D_p$와 projected oracle coefficient를 계산한다
([fusion_cate.py:168](../code/fusion_cate.py#L168)). 그러나 fresh draw가 20,000으로
고정되어 있다
([fusion_cate.py:176](../code/fusion_cate.py#L176)). v3는 moment별 MCSE를 확인하고
20k에서 40k, 80k, 최대 100k까지 늘리도록 한다
([roadmap:136](./implementation-experiment-roadmap.md#L136)).

또한 다음 oracle label을 결과와 그림에서 분리해야 한다.

- finite-grid true-risk oracle
- unregularized spline의 population sieve-moment oracle
- adaptive Monte Carlo feasible plug-in
- neural finite-grid oracle

### 4.6 Comparator fairness

**RCT-only 판정: `PASS`, exploratory.** 같은 basis와 score-based ridge 선택을 쓴다.

**Domain-indicator pooling 판정: `PASS`, 중요한 negative comparator.** Current synthetic
결과에서는 이 comparator가 FusionPPI보다 좋다. 이를 숨기지 않고 primary comparator로
유지한다
([monitoring report:242](./2026-09-12-implementation-monitoring-report.md#L242)).

**Experimental grounding 판정: `PASS`, heuristic label 필요.** neural에서 RCT-only보다
크게 나쁘고 spline에서는 거의 같다
([monitoring report:279](./2026-09-12-implementation-monitoring-report.md#L279)). 원 논문의
exact estimator라고 부르지 않고 mechanism-matched heuristic이라고 표시한다.

### 4.7 CATE-2 identification

**판정: `FAIL`, v3 source-law implementation.** Current 27-pattern construction은 그 coarse
law 안에서는 source $X$ law를 같게 만든다. 하지만 v3 exact tuple과 OOF freeze 계약을
만족하지 않는다.

**판정: `PASS`, exploratory reporting score.** spline은 두 outcomes와 두 losses에서
negative score gap을 보인다. neural은 mathematics와 reading의 방향이 다르다
([monitoring report:341](./2026-09-12-implementation-monitoring-report.md#L341)). 이 차이는
real outcome stress evidence다. source-law equality만으로 treatment contrast가 causal
CATE를 식별하지는 않는다
([handoff:846](./2026-09-09-implementation-handoff-JA.md#L846)).

### 4.8 CATE-3 transport

**판정: `PASS`, exploratory ratio sensitivity.** strong shift, $n_R=200$에서 true-ratio
selected/RCT median ratio는 spline $0.373$, neural $0.758$이다. unit ratio는 각각
$0.405$, $0.806$이다
([monitoring report:297](./2026-09-12-implementation-monitoring-report.md#L297)). 이것은
해당 cells에서 exact ratio가 unit보다 좋다는 결과다. “shift cost가 0이다” 또는 fitted
ratio theory가 confirm됐다는 뜻은 아니다.

### 4.9 Tail failures

Current headline median은 대형 위험 사건을 숨긴다. 다음은 current CSV를 read-only로
재집계한 **사후 audit diagnostic**이다. confirmatory success threshold가 아니다.

| file | diagnostic | result | unique seeds |
| --- | --- | ---: | --- |
| CATE-1 spline | selected/RCT risk ratio $>10$ | 3 learner rows, 2 seeds | `528146130`, `2897696311` |
| CATE-1 neural | selected/RCT risk ratio $>10$ | 4 learner rows, 4 seeds | `731449525`, `814266256`, `1737011207`, `3356062967` |
| CATE-1 spline | absolute selected risk $>100$ | 14 learner rows, 8 seeds | `313576480`, `490039293`, `528146130`, `607609814`, `684897006`, `2641424924`, `2897696311`, `3719665634` |
| CATE-1 neural | absolute selected risk $>100$ | 7 learner rows, 5 seeds | `38969307`, `1737011207`, `2030565177`, `2332066056`, `2392071668` |

두 diagnostic의 scale이 다르므로 서로 합치지 않는다. 이 사건은 실패 replication을
삭제할 근거가 아니다. seed별로 basis condition number, propensity extremes, selected
coefficient, candidate prediction quantiles, RCT-only risk를 재현하고 failure mechanism을
분류해야 한다.

## 5. Root causes

### 5.1 DGP가 level confounding만 만든다

Current hidden confounder는 outcome에 additively 들어가고 $X$와 독립적이다. OBS contrast
bias는 크지만 $X$에 따른 모양 변화가 작다. 예를 들어 Gaussian family 1,
confounding 2에서 mean bias는 $0.557$이고 leading-$X$ quintile spread는 $0.030$이다
([monitoring report:257](./2026-09-12-implementation-monitoring-report.md#L257)). source
indicator가 이 거의 상수인 차이를 흡수한다.

### 5.2 두 channel의 기회가 독립적으로 조절되지 않는다

현재 DGP는 $\lambda$ gain, $\omega$ gain, representation approximation, confounding shape를
동시에 바꾼다. 따라서 joint가 좋아져도 어느 channel이 왜 기여했는지 식별하기 어렵다.

### 5.3 평균적 성능과 tail safety를 분리하지 않는다

median risk ratio는 typical replication을 잘 보여주지만 rare instability를 숨긴다.
confirmatory design은 mean risk, median ratio, upper-tail ratio, catastrophic failure rate를
동시에 보고해야 한다.

### 5.4 v3와 rev4-style exploratory pipeline이 병존한다

현재 full CATE CSV는 v3의 four-role sampling, exact theorem path, adaptive oracle를 그대로
실행한 결과가 아니다. 이것을 code bug라고만 부르면 안 된다. 일부는 명시적인 protocol
deviation이고 일부는 추가 exploratory extension이다.

## 6. Manuscript-faithful method target

### 6.1 Fusion losses

DRF loss는

$$
\widehat L_{\mathrm{DRF}}(\zeta;\lambda,\omega,r)
=\mathbb P_R^{\mathrm{tune}}(Z_\lambda-\zeta)^2
+\omega\left[
\mathbb P_O^{\mathrm{tune}}\{r(\widehat g-\zeta)^2\}
-\mathbb P_R^{\mathrm{tune}}(\widehat g-\zeta)^2
\right].
$$

RF loss는 RCT R-loss term과 같은 transported squared-distance correction을 사용한다
([manuscript 4.tex:47](../manuscript/main/4.tex#L47)). Exact transport에서는 correction의
기대값이 0이므로 $\omega$는 target을 바꾸지 않고 empirical fitting variance를 바꾼다.

### 6.2 Sieve risk와 oracle channels

$b(X)\in\mathbb R^p$를 nuisance stage에서 고정한 basis라 하고
$\zeta_p=b^\top\beta_p$를 true CATE의 target-law projection이라 하자. DRF first-order
risk difference는

$$
n_R^{\mathrm{tune}}
\{\mathcal J_{\mathrm{DRF}}(\lambda,\omega)-\mathcal J_{\mathrm{DRF}}(0,0)\}
=2C_p\lambda-2D_p\omega+A_p\lambda^2+B_p\omega^2.
$$

따라서

$$
\lambda_p^\star=\Pi_{[0,1]}\left(-\frac{C_p}{A_p}\right),
\qquad
\omega_p^\star=\Pi_{\Omega}\left(\frac{D_p}{B_p}\right).
$$

정의와 gain 식은
[manuscript 4.tex:333](../manuscript/main/4.tex#L333)에 있다. RF는 일반적으로
$\lambda\omega$ cross term이 있고, $e(X)\equiv1/2$이면 그 term이 사라진다
([manuscript 4.tex:359](../manuscript/main/4.tex#L359)). 네 quadrant의 가장 깨끗한
mechanism check는 먼저 DRF와 constant trial propensity에서 수행한다. 이후 RF로
일반화한다.

## 7. Fair DGP suite

### 7.1 공통 구조

bounded core는 다음 구조를 사용한다.

$$
X\sim\operatorname{Unif}([-1,1]^d),
\qquad
Y(a)=m(X)+a\tau(X)+\eta(X,U)+\sigma(X)\epsilon,
$$

$$
\operatorname{logit}P_O(A=1\mid X,U)
=\operatorname{logit}e_R(X)+c_0U+c_1h(X)U.
$$

$U$는 learner가 보지 못하는 confounder다. $c_0U$는 거의 constant contrast bias를,
$c_1h(X)U$는 $X$-dependent bias를 만든다.

### 7.2 Four mechanism quadrants

quadrant는 이름이 아니라 independently estimated oracle moments로 확인한다.

| quadrant | moment condition | 필요한 결과 |
| --- | --- | --- |
| L0O0 | $C_p\approx0$, $D_p\approx0$ | 두 channel 모두 RCT fallback에 가까움 |
| L1O0 | $-C_p/A_p$가 interior, $D_p\approx0$ | $\lambda$만 유용함 |
| L0O1 | $C_p\approx0$, $D_p/B_p$가 interior | $\omega$만 유용함 |
| L1O1 | 두 oracle가 interior | joint gain의 기회가 있음 |

`approximately zero` tolerance는 outcome을 보기 전에 oracle-moment MCSE의 배수로
정한다. 예를 들어 zero interval이 0을 포함하는지로 판정하며 임의 absolute cutoff를
사후 선택하지 않는다.

### 7.3 Pure fixed-coordinate ablations

각 replication에서 다음 비교를 함께 계산한다.

1. **Pure $\lambda$ effect:** 같은 $\omega=0$에서
   $(\lambda^\star,0)$ 대 $(0,0)$
2. **Pure $\omega$ effect:** 같은 $\lambda=0$에서
   $(0,\omega^\star)$ 대 $(0,0)$
3. **Incremental $\omega$ effect:** 같은 fixed $\lambda$에서
   $(\lambda,\omega^\star)$ 대 $(\lambda,0)$
4. **Incremental $\lambda$ effect:** 같은 fixed $\omega$에서
   $(\lambda^\star,\omega)$ 대 $(0,\omega)$
5. **Selection effect:** selected coefficient 대 같은 candidate grid의 true-risk oracle

Oracle-coordinate ablation과 learned-coordinate ablation을 분리한다. Oracle result는
mechanism이 존재하는지 묻고, learned result는 selection procedure가 그 기회를 찾는지
묻는다.

### 7.4 Positive regimes

- L1O0, low nuisance bias와 informative $\Delta$
- L0O1, centered predictive $\widehat g$와 large $N_O/n_R$
- L1O1, 두 moment ratio가 모두 interior
- smooth $\tau$가 spline span에 가까운 regime
- $X$-dependent confounding을 domain indicator가 흡수하지 못하지만 FusionPPI의
  RCT correction과 selection이 통제할 수 있는 regime

### 7.5 Neutral controls

- L0O0
- $\widehat g$가 constant 또는 target-law에서 $\tau$와 orthogonal
- OBS outcome regression과 RCT regression의 prediction quality가 같아 $\Delta$가
  유용하지 않은 regime
- no confounding
- common marginal, $r_0=1$
- constant observational contrast bias, domain indicator에 유리한 regime

### 7.6 Adversarial and failure regimes

- $h(X)$가 sign-changing하는 confounding
- spline span보다 높은 frequency의 $h(X)$ 또는 $\tau(X)$
- adversarial $\widehat g$가 $\tau$와 음의 상관을 갖는 regime
- weak trial overlap
- heavy-tailed residuals
- small $n_R$와 ill-conditioned basis
- density-ratio misspecification
- ESS가 낮은 covariate shift
- ratio clipping이 target drift를 만드는 regime
- OBS pure noise 또는 severe outcome-model misspecification

Failure regime는 proposed method만 공격하지 않는다. 같은 data에서 RCT-only, single
channels, joint, grid oracle, domain-indicator pooling, experimental grounding을 모두
평가한다.

### 7.7 Transport panel

Gaussian mean shift는 analytic common support와 unbounded ratio를 사용한다. mild와
strong은 $d^2=\log2,\log5$로 고정하여

$$
\mathbb E_O(r_0^2)=e^{d^2},
\qquad
\operatorname{ESS\ fraction}=e^{-d^2}
$$

가 각각 $1/2,1/5$가 되게 한다
([roadmap:139](./implementation-experiment-roadmap.md#L139)). Five routes는 oracle,
classifier, BAL-X, BAL-X+$g$, unit이다. BAL-X+$g$를 표와 그림에서 생략하지 않는다.

## 8. Implementation PRD

### 8.1 P0, confirmatory validity

#### P0-1 Four-role CATE engine

- nuisance, tuning, selection, reporting을 각 source에서 독립 생성
- candidate는 tuning에서 fit, selection에서 선택, reporting 전에 freeze
- reporting outcome perturbation이 selection과 coefficient를 바꾸지 않아야 함
- role별 canonical RNG namespace와 seed fingerprint 기록

Acceptance:

- role row IDs disjoint 또는 independent source draw fingerprint가 다름
- 같은 seed 재실행은 scientific artifacts byte-identical
- reporting-only perturbation에서 selected index와 coefficient 변화 0
- selection-only perturbation에서 reporting raw sample 변화 0

#### P0-2 Theorem path separation

- Theorem 6: bounded clipped SCM-B, fixed $B_0=9$, $\bar r_0=1$
- Theorem 7/Corollary 7.1: spline, $\rho=0$, unclipped, $p<n_R^{tune}$,
  condition number $\le10^8$
- neural: empirical finite-grid extension label

Acceptance:

- theorem table에 neural과 clipped spline row가 0개
- theorem-7 row의 $\rho=0$과 unclipped flag가 모두 true
- condition failure는 explicit failure이며 redraw하지 않음
- bounded radius와 `bound_vacuity_ratio`를 함께 기록

#### P0-3 Adaptive oracle

- $A_p,B_p,C_p,D_p$를 20k에서 시작
- moment별 batch MCSE로 40k, 80k, 최대 100k까지 확장
- near-zero covariance에는 absolute MCSE 사용
- cap에서 불안정하면 `MC_INCONCLUSIVE`

Acceptance:

- 각 moment에 draws, estimate, MCSE, tolerance, status 존재
- oracle ordering QA는 oracle가 정의된 exact-ratio theorem cells에만 적용
- 같은 Monte Carlo sample로 만든 self-referential check 금지

#### P0-4 STAR v3 source law

- exact learner tuple은 gender, ethnicity3, free lunch, school type, school ID,
  corresponding missingness indicators
- exact cell $e_x\in[0.15,0.85]$
- five-fold deterministic OOF pooled outcome regression
- stable row-ID tie-break로 within-arm residual rank
- outcome·$\alpha$별 law를 replication 전에 한 번 freeze
- source hash, ordered row-ID hash, fold hash, probability table hash, propensity table
  hash 기록

Acceptance:

- $\max_x|P_O^X(x)-P_R^X(x)|\le10^{-12}$ at source-law level
- $\max_x|e(x)-P_R(A=1\mid X=x)|\le10^{-12}$
- OOF rank와 source weight가 learner feature에 없음
- nuisance, tune, select, report RNG streams가 독립
- source-law freeze 실패 시 hard G0 failure

#### P0-5 Tail failure ledger

- every scheduled seed에 result 또는 explicit failure
- risk와 risk ratio의 50%, 90%, 95%, 99%, maximum 보고
- prespecified catastrophic threshold와 numerical-failure threshold를 구분
- 현재 알려진 seed를 deterministic regression fixtures로 보존

Acceptance:

- scheduled seed count = successful seed count + failed seed count
- duplicate 또는 redrawn seed 0
- 현재 표의 대형 위험 seed를 모두 재현하거나 code change에 따른 차이를 설명

### 8.2 P1, scientific endpoints

#### P1-1 Paired risk inference

- cell 내부에서 method difference를 replication-paired로 계산
- family aggregation은 family-stratified paired bootstrap
- point estimate와 independent replication count 기록

#### P1-2 Variance and squared bias

fixed target-law anchors에서

$$
\overline R_K
=\frac{K-1}{K}\widehat V_K+\widehat B_K^2
$$

를 사용한다. $K$는 independent replications 수다. integrated anchor Monte Carlo error를
별도로 평가한다.

Acceptance:

- decomposition gap $\le10^{-10}$ on deterministic fixture
- stochastic result에는 anchor size와 integration MCSE 존재
- variance, squared bias, total risk의 scale과 aggregation이 동일

#### P1-3 Score calibration and regret

- candidate-level score difference 대 true-risk difference
- calibration slope/intercept와 rank correlation
- selected score optimism
- selected-to-grid-oracle regret
- best-single-channel과 RCT-only 대비 excess risk

Acceptance:

- selection sample과 reporting sample 결과를 별도 column으로 저장
- true-risk column이 STAR row에는 `N/A`
- regret가 negative이면 oracle label 또는 candidate set 오류로 hard fail

#### P1-4 Coefficient mechanism

- selected, grid oracle, sieve-moment oracle, feasible plug-in coefficient를 분리
- $\lambda$, $\omega$, $\rho$ boundary frequency
- denominator instability와 fallback frequency
- fixed-coordinate ablation 결과

Acceptance:

- every selected candidate가 declared grid에 포함
- deterministic tie-break가 grid order와 일치
- reporting perturbation으로 coefficient가 움직이지 않음

#### P1-5 Ratio diagnostics

- $\mathbb E_O\widehat r$, ESS, $\mathbb E_O\widehat r^2$
- log-ratio quantiles
- BAL feature residual
- normalization error
- target drift diagnostic
- raw optimizer status와 protocol balance status 분리

Acceptance:

- BAL-X와 BAL-X+$g$ residual $\le10^{-6}$일 때만 protocol-converged
- non-BAL balance residual은 `N/A`
- normalization error는 diagnostic이고 arbitrary hard gate가 아님
- all five routes가 summary와 plot에 존재

### 8.3 P2, extensions

- pointwise CATE interval와 simultaneous bands
- pointwise coverage와 interval width
- policy value와 policy regret
- survival or censoring extension for WHI
- estimated-ratio CATE theorem과 cross-fitted reuse theory

P2는 현재 point-estimation manuscript를 검증하는 필수 구현이 아니다. P0와 P1이
완료되기 전에 P2를 이유로 confirmatory run을 지연하지 않는다.

## 9. Agile execution M0 to M5

### M0. Evidence freeze and protocol lock

산출물:

- existing CATE artifact inventory와 SHA-256 manifest
- code, manuscript, handoff hash
- exact schedule와 seed ledger
- exploratory/confirmatory label map

Go:

- current artifacts를 덮어쓰지 않는 atomic stage promotion이 검증됨
- code 또는 DGP 변경 시 prefix incompatibility를 감지함

Stop:

- source artifact 또는 seed provenance를 복원할 수 없음

### M1. Deterministic mechanism micro fixtures

L0O0, L1O0, L0O1, L1O1 각 하나의 fixed-array fixture를 만든다. 이는 statistical
replication이 아니다.

Go:

- independent formula와 production $A_p,B_p,C_p,D_p$ 차이 $\le10^{-10}$
- oracle coordinate와 fixed-coordinate risk ordering이 fixture 설계와 일치
- deliberate sign, role, ratio mutation을 test가 검출

Stop and `REDESIGN`:

- identity, role, target, normal equation, oracle sign이 맞지 않음

### M2. One-rep end-to-end smoke

최소 path:

- four quadrant spline DRF
- one RF cross-term fixture
- bounded Theorem 6
- unclipped $\rho=0$ Theorem 7
- neural empirical path
- five transport routes
- STAR mathematics $\alpha=0.8,n_R=200$ source law and reporting path

Go:

- P0 hard assertions 전부 통과
- every path가 result 또는 explicit expected failure 생성

Stop and `REDESIGN`:

- source law, role independence, propensity, formula, provenance failure

Directional performance failure는 engineering failure가 아니다.

### M3. Ten-rep pilot

| panel | cells | reps/cell | independent reps |
| --- | ---: | ---: | ---: |
| core quadrants × 2 confounding shapes × 3 $n_R$ | 24 | 10 | 240 |
| O1/L1O1 × 2 shifts × 3 $n_R$ | 12 | 10 | 120 |
| six negative/failure controls | 6 | 10 | 60 |
| 합계 | 42 |  | 420 |

목적:

- finite output와 failure mechanism
- direction과 effect scale
- learner·sieve·route runtime
- confirmatory power와 MCSE planning

10회로 superiority, coverage, no-harm을 주장하지 않는다. Confirmatory outcome을 열기
전에 primary contrast, interval convention, multiplicity adjustment, minimum material
effect $\delta$, noninferiority margin을 freeze한다. 기존 gate의 5%는 protocol choice이지
통계 법칙이 아니다.

### M4. Decision-directed escalation

필요한 inconclusive mechanism cells만 50회까지 늘린다. 모든 cell을 자동으로 채우지
않는다.

- validity failure: `REDESIGN`
- valid implementation이지만 channel opportunity 없음: `REFRAME`
- valid하고 방향이 맞지만 uncertainty가 큼: `INCONCLUSIVE`
- prespecified scientific criteria를 만족: `GO_CONFIRMATORY`

사후에 DGP 파라미터를 바꾼 결과는 새 exploratory version으로 분리한다. 실패 seed를
교체하지 않는다.

### M5. Confirmatory simulation

| panel | cells | reps/cell | independent reps |
| --- | ---: | ---: | ---: |
| four quadrants × 2 confounding shapes × 3 $n_R$ | 24 | 200 | 4,800 |
| O1/L1O1 × 2 shift levels × 3 $n_R$ | 12 | 200 | 2,400 |
| six negative/failure controls | 6 | 100 | 600 |
| 합계 | 42 |  | 7,800 |

7,800은 independent DGP replications다. learner, sieve, route, candidate, metric row를
곱한 model-fit count가 아니다. M3 prefix는 protocol, code path, seed가 완전히 같을 때만
재사용한다. 불명확하면 기존 pilot을 보존하고 confirmatory IDs를 새로 실행한다.

Confirmatory go criteria:

- L1O0에서 fixed $\omega=0$의 $\lambda$ contrast가 prespecified criterion 충족
- L0O1에서 fixed $\lambda=0$의 $\omega$ contrast가 criterion 충족
- L1O1에서 joint가 best single channel보다 개선
- L0O0에서 prespecified noninferiority와 fallback behavior 충족
- varying-bias cells에서 domain-indicator comparator와의 primary paired comparison 보고
- tail failure rate가 prespecified safety criterion 충족
- 모든 primary interval과 multiplicity rule이 사전 동결됨

Joint가 RCT-only를 이기지 못한 것 자체는 code failure가 아니다.

## 10. Endpoints

### 10.1 Primary synthetic endpoint

cell별 paired mean risk difference:

$$
\Delta_R
=\mathbb E_{\mathrm{rep}}
[\mathcal R(\widehat\zeta_{\mathrm{Fusion}})
-\mathcal R(\widehat\zeta_{\mathrm{RCT}})].
$$

음수면 FusionPPI의 평균 true risk가 낮다. interval은 replication-paired이고 family
aggregation은 family-stratified다.

### 10.2 Primary mechanism contrasts

- same-$\omega$ $\lambda$ effect
- same-$\lambda$ $\omega$ effect
- joint 대 best feasible single channel
- selected 대 grid oracle regret
- FusionPPI 대 domain-indicator pooling

### 10.3 Secondary performance endpoints

- absolute integrated risk
- per-rep risk ratio와 median
- integrated prediction variance
- integrated squared bias
- total anchor risk와 decomposition gap
- 90%, 95%, 99% risk와 risk-ratio quantiles
- catastrophic and numerical failure rates
- selection regret
- selected coefficient distribution and stability
- score-to-risk calibration
- ratio ESS, balance, normalization, drift
- runtime and memory

Mean risk difference가 primary다. Ratio는 RCT risk가 매우 작을 때 불안정하므로 paired
difference와 함께 해석한다.

### 10.4 Real-data endpoints

STAR, NSW, WHI에는 pointwise true $\tau(X)$가 없다. Primary endpoint는 independent RCT
reporting-score difference다. Secondary endpoints는 descriptive RCT ATE benchmark,
prediction stability, policy value가 가능한 경우의 experimental value다. 이를 true CATE
risk 또는 unbiased performance minimum이라고 부르지 않는다.

### 10.5 Pointwise coverage

향후 interval extension에서만

$$
\operatorname{Coverage}(x)
=P\{\tau(x)\in[\widehat L(x),\widehat U(x)]\}
$$

를 평가한다. average pointwise coverage, worst-subgroup coverage, simultaneous coverage를
분리한다. 현재 integrated risk와 variance decomposition만으로 coverage를 주장하지 않는다.

## 11. Real-world plan

### 11.1 Repaired STAR primary stress test

#### Estimand

frozen eligible kindergarten cohort에서 small class 대 pooled regular and regular-with-aide
control의 mathematics primary, reading secondary CATE다.

#### Source construction

- $P_R$: frozen eligible cohort의 uniform empirical law
- $P_O$: 먼저 exact $X$를 $P_R^X$에서 뽑고 $A$를 exact-cell propensity에서 뽑은 뒤,
  같은 $(X,A)$ cell에서 OOF outcome-rank weight로 row를 뽑음
- primary ratio: $r_0=1$ at source-law level
- secondary covariate-shift STAR panel은 별도 prospective protocol

#### Roles

- RCT nuisance and tuning: $n_R\in\{200,400\}$ each
- OBS nuisance and tuning: $N_O=5000$ each
- selection and reporting: source별 1,000 each

이는 v3 exact design이다
([handoff:904](./2026-09-09-implementation-handoff-JA.md#L904)).

#### Validation

- independent reporting-score difference
- eligible-cohort descriptive ATE
- score calibration across replicated source draws
- true CATE risk claim 금지
- source-law equality가 causal identification을 증명한다는 claim 금지

#### Reproducibility and licensing

- official provenance와 data dictionary 확인
- source URL, exact file SHA-256, row count, schema freeze
- 현재 raw GitHub CSV의 redistribution license를 확인하기 전에는 archive 또는 배포 금지

### 11.2 NSW with CPS and PSID, supporting benchmark only

#### Estimand

NSW experimental target population에서 job training이 1978 earnings에 미치는 heterogeneous
effect다.

#### Construction

- RCT: NSW experimental treated and controls
- OBS comparison: CPS 또는 PSID controls
- target ratio: experimental target law 대 comparison covariate law
- treated unit 재사용을 피하도록 역할별 partition 또는 prospective cross-fitting 사용

#### Validation

- independent experimental reporting score
- experimental ATE benchmark
- prespecified policy value가 가능하면 secondary
- known pointwise CATE truth와 oracle ratio를 꾸미지 않음

#### Position

표본이 작고 observational controls가 역사적 benchmark이므로 primary confirmatory real-data
claim이 아니라 supporting stress test다. CPS와 PSID를 모두 보고한다.

#### Reproducibility

- original NBER or package version, URL, SHA-256, row count
- exact variable and treatment/outcome crosswalk
- redistribution terms 확인

### 11.3 WHI CT plus OS, strongest natural-pair feasibility

#### Candidate estimand

WHI Clinical Trial eligibility에 맞춘 target population에서 hormone therapy의 5-year binary
risk difference CATE 또는 restricted mean survival time CATE다. 둘 중 하나를 access 후
outcome inspection 전에 고정한다.

#### Construction

- RCT: WHI Clinical Trial
- OBS: WHI Observational Study
- baseline time, eligibility, treatment definition, outcome definition을 semantic crosswalk로
  일치시킴
- ratio: $dP_{CT}^X/dP_{OS}^X$
- participant 또는 site cluster를 분리한 nuisance, tune, select, report roles

#### Validation

- held-out CT reporting score
- marginal CT ATE benchmark
- survival outcome이면 censoring nuisance와 estimand-specific loss를 먼저 구현

#### Feasibility gate

- participant-level data access와 DUA
- CT/OS common-variable crosswalk
- treatment initiation and follow-up alignment
- event count와 overlap
- censoring support
- permitted artifact disclosure

WHI는 자연적으로 병렬인 RCT와 OBS라는 점에서 가장 강한 real-world 후보지만, access와
survival extension이 끝나기 전에는 실행 예산을 약속하지 않는다.

## 12. Safe and forbidden claims

### 12.1 Safe now

- Current exploratory synthetic cells에서 selected FusionPPI가 fair RCT-only보다 낮은
  median risk를 보인 경우가 많았다.
- Current gain은 반복 fitting의 prediction variance 감소와 관련되어 있다
  ([monitoring report:232](./2026-09-12-implementation-monitoring-report.md#L232)).
- Current synthetic DGP에서는 domain-indicator pooling이 더 좋다.
- Strong Gaussian shift cells에서 exact ratio route가 unit route보다 낮은 median risk를
  보였다.
- STAR spline score는 개선됐지만 neural 결과는 outcome에 따라 다르다.

### 12.2 Forbidden until confirmatory completion

- “FusionPPI는 observational confounding에서 pooling보다 일반적으로 우월하다.”
- “Estimated coefficients는 finite sample no-harm을 보장한다.”
- “Known ratio면 covariate shift 비용이 0이다.”
- “STAR에서 true CATE risk를 줄였다.”
- “Source-law equality가 STAR causal identification을 증명한다.”
- “Neural result가 Theorem 7 또는 Corollary 7.1을 검증한다.”
- “CATE-3가 current common-marginal theorem을 검증한다.”
- “Median improvement가 catastrophic tail safety를 보장한다.”
- “Grid oracle 또는 moment oracle가 실행 가능한 estimator다.”
- “Low integrated risk가 pointwise interval coverage를 보장한다.”

## 13. Plot PRD

모든 plot은 독립 replication count, aggregation estimand, uncertainty unit을 caption에
쓴다. route-row count를 replication count라고 부르지 않는다.

### Figure C1. Channel opportunity map

- x-axis: $\lambda_p^\star$
- y-axis: $\omega_p^\star$
- color: quadrant and confounding shape
- point interval: oracle MCSE propagated interval
- annotation: $A_p,B_p,C_p,D_p$, sample sizes
- purpose: DGP가 실제로 L0O0/L1O0/L0O1/L1O1을 만들었는지 확인

### Figure C2. Primary paired risk differences

- panels: DRF and RF, spline and neural
- y-axis: paired mean $\Delta_R$
- interval: family-stratified paired bootstrap 95%
- rows: four quadrants and negative controls
- comparators: RCT-only, $\lambda$-only, $\omega$-only, joint, domain indicator
- zero reference line

### Figure C3. Variance and squared-bias decomposition

- stacked components: integrated prediction variance and squared bias
- normalize within sieve by matched RCT total, so RCT bar is 1
- annotate raw total risk
- show selected and grid oracle separately
- include decomposition gap in companion table, not rounded away

### Figure C4. Coefficient recovery and selection regret

- panel A: selected versus moment-oracle $\lambda$
- panel B: selected versus moment-oracle $\omega$
- panel C: boundary and fallback rates
- panel D: selection regret distribution
- distinguish grid oracle and moment oracle by symbol and caption

### Figure C5. Transport routes

- routes: oracle, classifier, BAL-X, BAL-X+$g$, unit
- no route omission
- panels: mild/strong, spline/neural, DRF/RF
- primary y-axis: paired mean risk difference versus RCT-only
- secondary table: ESS, balance residual, normalization error, drift
- Gaussian ratio is unbounded, so no sample-maximum global bound

### Figure C6. Tail and failure diagnostics

- empirical CDF or complementary CDF of risk ratio
- 90%, 95%, 99%, maximum markers
- catastrophic threshold prespecified before confirmatory outcomes
- separate numerical failure, protocol failure, scientific negative transfer
- label current known seeds in exploratory appendix only

### Figure C7. Real-data reporting score

- STAR mathematics and reading separate
- NSW CPS and PSID separate
- WHI only after feasibility gate
- y-axis: independent reporting score selected minus basis-matched RCT-only
- zero reference line
- never label as true CATE risk

## 14. Table PRD

### Table C1. Design and provenance

cell ID, DGP version, sample role sizes, learner, sieve, ratio, replication IDs, code hash,
protocol hash, source hash, failures.

### Table C2. Primary risk results

paired mean difference, 95% interval, absolute risks, median ratio, independent replication
count, multiplicity status.

### Table C3. Mechanism results

$A_p,B_p,C_p,D_p$, MCSE, oracle coefficients, selected coefficients, pure fixed-coordinate
contrasts, fallback rate.

### Table C4. Variance and bias

integrated prediction variance, squared bias, total risk, anchor size, integration MCSE,
decomposition gap.

### Table C5. Calibration and regret

score-risk slope, intercept, rank correlation, selection optimism, mean and tail regret.

### Table C6. Ratio and overlap

five routes, $E_O\widehat r$, $E_O\widehat r^2$, ESS, log-ratio quantiles, balance residual,
normalization error, drift, optimizer status.

### Table C7. Tail and failures

failure category, count, denominator, rate, seed IDs, resolved or unresolved mechanism.

### Table C8. Real-data results

dataset, estimand, source construction, reporting-score difference, descriptive ATE, overlap,
license/access status, allowed claim.

## 15. Acceptance-test matrix

| ID | requirement | fixture or data | observable pass condition | failure class |
| --- | --- | --- | --- | --- |
| T01 | pseudo-outcome validity | fixed bounded arrays | independent conditional-mean identity within $10^{-10}$ | hard |
| T02 | DRF normal equation | fixed arrays | direct objective minimizer and solve agree within $10^{-10}$ | hard |
| T03 | RF normal equation | fixed arrays | direct objective and solve agree within $10^{-10}$ | hard |
| T04 | four role separation | generated IDs/fingerprints | nuisance, tune, select, report independent | hard |
| T05 | reporting freeze | perturb reporting only | candidate index and coefficients unchanged | hard |
| T06 | selection sensitivity | perturb selection only | reporting sample unchanged | hard |
| T07 | canonical seeds | namespace fixture | stable, distinct 128-bit-derived streams | hard |
| T08 | exact common marginal | source-law table | maximum probability difference $\le10^{-12}$ | hard |
| T09 | exact propensity | STAR exact cells | maximum identity error $\le10^{-12}$ | hard |
| T10 | OOF residual | STAR folds | every row predicted by fold excluding that row | hard |
| T11 | no design leakage | STAR feature schema | rank, residual, weight absent from learner $X$ | hard |
| T12 | theorem scope | result schema | only eligible unclipped $\rho=0$ spline in theorem table | hard |
| T13 | bounded theorem | raw SCM-B components | outcome, score, prediction within prespecified bounds | hard |
| T14 | radius vacuity | theorem rows | radius and vacuity ratio both finite and reported | report |
| T15 | sieve oracle moments | fixed arrays | production and independent $A_p$ to $D_p$ agree $\le10^{-10}$ | hard |
| T16 | adaptive oracle | batch fixture | draws and MCSE stopping rule exactly followed | hard |
| T17 | grid oracle ordering | exact-ratio synthetic | oracle risk no larger than every candidate up to $10^{-12}$ | hard |
| T18 | selection regret | synthetic | regret $\ge-10^{-12}$ | hard |
| T19 | risk decomposition | fixed anchor matrix | gap $\le10^{-10}$ | hard |
| T20 | score identity | discrete support | exact score-risk difference identity within $10^{-10}$ | hard |
| T21 | score calibration | synthetic replication set | all candidate scores and risks paired by ID | hard |
| T22 | BAL-X convergence | feasible fixture | maximum residual $\le10^{-6}$ | hard |
| T23 | BAL-X+$g$ content | feature mutation | route differs from BAL-X and balances declared feature | hard |
| T24 | five-route completeness | schedule validator | exactly oracle/classifier/BAL-X/BAL-X+$g$/unit | hard |
| T25 | Gaussian moments | analytic fixture | $E_Or_0=1$, $E_Or_0^2=e^{d^2}$, ESS $=e^{-d^2}$ | hard |
| T26 | no global ratio bound | language/schema scan | Gaussian radius N/A, no sample-max bound claim | hard |
| T27 | failure accounting | schedule manifest | scheduled = success + explicit failure | hard |
| T28 | no redraw | seed ledger | duplicate/replacement seed count 0 | hard |
| T29 | current tail seeds | regression subset | known exploratory events reproduced or explained | diagnostic |
| T30 | bootstrap pairing | fixed sample | manual paired calculation and implementation agree | hard |
| T31 | independent rep count | summary schema | count uses family-cell-rep IDs, not rows | hard |
| T32 | runtime provenance | pilot | role/learner/sieve/route runtime and code hash present | hard |
| T33 | STAR metric label | real-data rows | true-risk fields N/A, reporting-score label present | hard |
| T34 | CATE-3 theory label | result/report scan | `exploratory estimated-ratio extension` present | hard |
| T35 | pointwise coverage label | inference extension only | no coverage claim without constructed intervals | hard |

90% 이상의 PRD 요구가 위와 같이 자동 test 또는 machine-readable validator로 확인돼야
한다. “reviewer가 읽어 보면 맞다”만으로 acceptance를 대신하지 않는다.

## 16. Compute budgeting

runtime은 model fit 수에 따라 계산한다. 독립 DGP 반복 수만으로 벽시계 시간을 단순
외삽하지 않는다.

$$
T_{\mathrm{single}}
=\sum_{c\in\mathcal C}
\sum_{\ell\in\mathcal L}
\sum_{s\in\mathcal S}
\sum_{r\in\mathcal Q}
n_{c}\,\widetilde T(c,\ell,s,r),
$$

여기서 $\mathcal C$는 design cells, $\mathcal L$은 DRF/RF learners,
$\mathcal S$는 sieves, $\mathcal Q$는 ratio routes, $\widetilde T$는 M3에서 측정한 median
runtime이다.

M3 후 다음을 보고한다.

- learner·sieve·route별 median, 90%, maximum runtime
- peak resident memory
- oracle Monte Carlo share
- 1, 2, 4, 8 workers의 measured scaling
- failure retry가 아니라 failed-seed recording에 든 시간
- expected M5 wall time와 20% operational reserve

측정 전에는 8-core 시간이나 최종 비용을 약속하지 않는다. Neural이 병목이면 epoch를
조용히 줄이지 않는다. prospective protocol revision 또는 explicit failure로 처리한다.

## 17. Exact implementation order

1. M0 artifact와 provenance freeze
2. M1 four-quadrant fixed-array tests
3. P0-1 four-role engine
4. P0-2 theorem path separation
5. P0-3 adaptive oracle
6. P0-4 STAR source law
7. P0-5 schedule, failure, tail ledger
8. P1 endpoint and calibration implementation
9. M2 one-rep smoke
10. independent code, statistics, and manuscript-label review
11. M3 10-rep pilot
12. primary threshold, interval, multiplicity, runtime budget freeze
13. M4 decision-directed escalation
14. exactly one terminal decision record
15. `GO_CONFIRMATORY`일 때만 M5
16. M5 결과가 고정된 뒤 plots와 tables 생성
17. STAR repaired run
18. NSW/CPS and NSW/PSID supporting runs
19. WHI feasibility가 통과한 경우에만 WHI implementation plan
20. 별도 승인 후 manuscript update, mirror sync, commit, push

## 18. Final decision rule

### Engineering decision

- formula, role, source law, theorem label, provenance failure: `REDESIGN`
- failed seed를 삭제하거나 바꾼 경우: `REDESIGN`

### Scientific decision

- valid하고 사전 지정된 channel opportunity와 safety criteria를 충족:
  `GO_CONFIRMATORY`
- valid하지만 relevant cell이 cap 아래이고 uncertainty가 큼: `INCONCLUSIVE`
- valid하지만 mechanism opportunity 또는 broad claim이 지지되지 않음: `REFRAME`

세 scientific state와 engineering failure를 exhaustive하게 기록한다. Joint가 RCT-only보다
나쁜 scientific outcome은 자동 code failure가 아니다. 반대로 median이 좋다는 이유로
formula나 honesty failure를 무시하지 않는다.

## 19. Completion boundary

이 계획의 최소한의 완전한 목표는 FusionPPI가 도움을 받을 수 있는 조건, 도움을 받지
못하는 조건, 위험해지는 조건을 각각 공정하게 보여주는 것이다. “좋은 그림을 만드는 것”이
목표가 아니다. Positive, neutral, adversarial evidence를 같은 protocol에서 보존해야 한다.

현재 다음 단계는 M0과 M1의 구현계획을 별도 승인받는 것이다. 이 문서 자체는 code나
실험을 실행하지 않는다.

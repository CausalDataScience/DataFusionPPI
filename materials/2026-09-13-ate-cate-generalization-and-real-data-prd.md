# FusionPPI ATE-CATE 일반화 및 실제자료 실험 PRD

**상태:** 실행 전 설계 계약

**작성일:** 2026-09-13

**대상 독자:** 연구 책임자, 구현자, 결과 분석자, 원고 작성자

**통계 사양 authority:** `2026-09-09-implementation-handoff-JA.md`와 현재 원고

**이 문서의 역할:** 새 일반화 실험을 실행하기 위한 운영 계약

이 문서는 기존 통계 공식을 대체하지 않는다. 원고와 handoff가 estimand와
추정량의 통계적 authority이고, 이 문서는 데이터 생성, 표본 역할, 실행 순서,
산출물, 판정, 재개 규칙을 고정한다. 구현과 이 문서가 충돌하면 실행을 멈추고
충돌을 먼저 해결한다.

---

## 1. Output Contract

### 1.1 Deliverable

동일한 simulation bundle에서 ATE와 CATE를 함께 평가하는 다음 결과를 만든다.

1. 18개 synthetic primary cell의 독립 confirmatory 결과
2. STAR `mathk`에 대한 200개 conditional resplit 결과
3. ATE와 CATE를 분리해서 해석하는 그림 3개와 표 2개
4. 모든 실행을 재현하고 감사할 수 있는 code, protocol, input, seed, failure,
   environment, raw, summary, plot provenance

### 1.2 Consumer

연구 책임자는 이 결과로 FusionPPI가 특정 SCM이나 표본 크기에만 의존하는지,
두 fusion channel이 서로 다른 실패 조건에서 어떻게 작동하는지, 실제 outcome을
사용한 STAR에서 어떤 수준까지 안전하게 주장할 수 있는지를 판단한다. 원고
작성자는 사전 고정된 claim map을 넘어서는 문장을 쓰지 않는다.

### 1.3 Acceptance Criteria

완료는 유리한 결과가 나왔다는 뜻이 아니다. 다음 세 조건이 모두 충족되어야
한다.

1. 정확히 18개 synthetic cell을 각각 200회 실행하고 STAR conditional bundle을
   200회 실행한다. 운영 장부상 main source-data bundle은
   $3{,}600+200=3{,}800$개지만, synthetic과 STAR는 서로 다른 확률실험이므로
   추론에서 합치지 않는다.
2. ATE와 CATE가 각 replication에서 같은 원자료 bundle을 사용하되, nuisance,
   tuning, selection, evaluation 역할은 서로 독립이다.
3. 모든 예정 seed가 성공 결과 또는 명시적 failure record로 남고, raw 결과에서
   표와 그림이 재생성된다.

### 1.4 Evidence

완료 근거는 설명문이 아니라 immutable stage directory의 manifest, SHA256,
seed ledger, failure ledger, raw sufficient statistics, verification report,
summary table과 figure reproduction check이다.

### 1.5 Boundaries

이 PRD는 다음을 허용하지 않는다.

- 현재 원고의 theorem 또는 estimand 변경
- pilot 결과를 본 뒤 cell, coefficient grid, outcome perturbation 또는 threshold 변경
- 실패 replication 삭제 또는 새 seed로 교체
- 기존 MM-main 결과와 새 confirmatory 결과의 inferential pooling
- STAR를 causal ground truth 실험으로 표현
- 선언되지 않은 outcome, prediction 또는 density-ratio clipping
- WHI 자료 다운로드, 접근 신청 또는 분석
- dependency 추가, mirror 동기화, Git stage, commit 또는 push

---

## 2. Authority와 현재 상태

### 2.1 Authority 순서

1. [현재 원고의 표본·estimand 정의](../manuscript/main/2.tex#L7-L70)
2. [ATE 추정량과 분산식](../manuscript/main/3.tex#L14-L104)
3. [CATE 고정 basis와 weighted loss](../manuscript/main/4.tex#L230-L360)
4. [implementation handoff v3](2026-09-09-implementation-handoff-JA.md)
5. 이 실행 PRD
6. 실행 당시의 fingerprint-qualified protocol과 code manifest

상위 authority와 하위 문서가 충돌하면 상위 authority가 우선한다. 다만 실행
전에 충돌을 기록하고 해소해야 하며, 구현자가 임의로 해석해서 진행하지 않는다.

### 2.2 현재 존재하는 근거

기존 manuscript-minimal 실험은 SCM1, 역할별 $n_R=100$에서 shared와 shifted
두 cell을 각각 200회 실행했다. 결과에는 400개 cell record와 0개 failure가
기록되어 있다. 구체적 수치와 검증은
[MM-main RECORD](confirm/MM-main-66e93079cbd7/RECORD.md)에 있다.

이 결과는 새 실험의 설계 동기와 외부 일관성 비교에는 사용한다. 그러나 새
18-cell primary estimate, confidence interval 또는 multiplicity family에는 넣지
않는다. 기존 결과에 맞춰 새 설계를 만든 뒤 같은 결과를 다시 confirmatory
evidence로 세는 순환을 막기 위해서다.

기존 ATE 결과와 재현 manifest는
[ATE 결과 보고서](2026-09-12-ate-results-report.md)와
[ATE artifact manifest](ate-artifact-manifest.sha256)에 있다. 이것도 새 paired
replication과 합치지 않는다. 새 실험의 ATE는 새 CATE와 같은 bundle에서 다시
계산한다.

### 2.3 아직 존재하지 않는 것

다음은 planned implementation이며 현재 완료됐다고 주장하지 않는다.

- 18-cell 통합 synthetic schedule
- ATE와 CATE를 함께 생성하는 통합 role-bundle driver
- outcome-domain-mismatch DGP
- STAR ATE-CATE 통합 200-bundle driver
- 이 PRD의 그림 3개와 표 2개를 만드는 분석기

현재 `mm_design.py`, `mm_run.py`, `mm_analyze.py`는 재사용 후보인 검증된 구성요소다.
그 자체가 새 18-cell 실험의 완성된 실행기는 아니다.

---

## 3. 용어와 estimand

### 3.1 자료원

$R$은 randomized trial source, $O$는 observational source다. 한 관측치는

$$
V=(X,A,Y)
$$

로 쓴다. $X$는 공변량, $A\in\{0,1\}$은 처리, $Y$는 관측 outcome이다.

### 3.2 ATE와 CATE

RCT 모집단에서의 conditional average treatment effect, 즉 CATE는

$$
\tau(x)=\mathbb E_R\{Y(1)-Y(0)\mid X=x\}
$$

이다. Average treatment effect, 즉 ATE는

$$
\theta=\mathbb E_R\{\tau(X)\}
$$

이다. 원고의 정의는 [2.tex](../manuscript/main/2.tex#L25-L40)에 있다.

### 3.3 Density ratio

OBS 공변량 분포를 RCT 공변량 분포로 옮기는 exact density ratio는

$$
r_0(x)=\frac{dP_R^X}{dP_O^X}(x)
$$

이다. $P_R^X=P_O^X$이면 $r_0(x)=1$이다. Covariate shift에서는 support overlap이
필요하지만, overlap은 $r_0$가 알려졌거나 정확히 추정됐다는 뜻이 아니다.
이 구분은 [2.tex](../manuscript/main/2.tex#L42-L70)을 따른다.

### 3.4 Fusion coefficients

$\lambda$는 RCT와 OBS outcome regression을 RCT-valid score 안에서 혼합하는
계수다. $\omega$는 OBS prediction의 transported average와 RCT average의 차이를
control variate로 사용하는 계수다. ATE estimator는

$$
\widehat\theta_r(\lambda,\omega)
=
\mathbb P_R Z_\lambda
+\omega\left[
\mathbb P_O\{r(X)\widehat g(X)\}
-\mathbb P_R\widehat g(X)
\right]
$$

이다. Exact transport 아래의 표적과 분산식은
[3.tex](../manuscript/main/3.tex#L14-L85)을 따른다.

### 3.5 CATE risk

Synthetic experiment에서는 $\tau(X)$를 알고 있으므로 fitted CATE
$\widehat\tau(X)$의 integrated squared risk를 계산한다.

$$
R(\widehat\tau)
=
\mathbb E_R\left[
\{\widehat\tau(X)-\tau(X)\}^2
\right].
$$

이것은 ATE estimator의 MSE와 단위와 대상이 다르다. 두 수치를 한 축에서 직접
비교하지 않는다.

---

## 4. 연구 질문

이 실험은 다음 네 질문에 답한다.

1. SCM과 $n_R$이 바뀌어도 joint fusion의 ATE·CATE 성능 방향이 유지되는가?
2. Covariate shift에서 exact transport, estimated transport, 잘못된 transport가
   어떻게 다른가?
3. 공변량 분포는 같지만 OBS baseline outcome이 RCT와 다를 때 $\lambda$와
   $\omega$ channel은 어떻게 반응하는가?
4. 참 treatment-effect function을 모르는 STAR에서 held-out 성능과 coefficient
   behavior가 synthetic 결과와 양립하는가?

이 질문은 “모든 조건에서 fusion이 이긴다”는 결론을 전제로 하지 않는다.
성능 악화도 설계가 유효하면 과학적 결과다.

---

## 5. 통합 18-cell synthetic design

### 5.1 축

Primary synthetic rectangle은 다음 Cartesian product다.

$$
\mathrm{SCM}\in\{1,2,3\},\qquad
n_R\in\{100,400\},\qquad
q\in\{\text{shared},\text{covariate shift},\text{outcome mismatch}\}.
$$

따라서 cell 수는

$$
3\times2\times3=18
$$

개다. $n_R$은 각 RCT 역할에 들어가는 표본 수다.

### 5.2 정확한 cell 표

| cell_id | SCM | 역할별 $n_R$ | regime | exact ratio |
|---|---:|---:|---|---|
| S01 | 1 | 100 | shared | $r_0=1$ |
| S02 | 1 | 100 | covariate shift | analytic $r_0$ |
| S03 | 1 | 100 | outcome mismatch | $r_0=1$ |
| S04 | 1 | 400 | shared | $r_0=1$ |
| S05 | 1 | 400 | covariate shift | analytic $r_0$ |
| S06 | 1 | 400 | outcome mismatch | $r_0=1$ |
| S07 | 2 | 100 | shared | $r_0=1$ |
| S08 | 2 | 100 | covariate shift | analytic $r_0$ |
| S09 | 2 | 100 | outcome mismatch | $r_0=1$ |
| S10 | 2 | 400 | shared | $r_0=1$ |
| S11 | 2 | 400 | covariate shift | analytic $r_0$ |
| S12 | 2 | 400 | outcome mismatch | $r_0=1$ |
| S13 | 3 | 100 | shared | $r_0=1$ |
| S14 | 3 | 100 | covariate shift | analytic $r_0$ |
| S15 | 3 | 100 | outcome mismatch | $r_0=1$ |
| S16 | 3 | 400 | shared | $r_0=1$ |
| S17 | 3 | 400 | covariate shift | analytic $r_0$ |
| S18 | 3 | 400 | outcome mismatch | $r_0=1$ |

모든 cell을 새 confirmatory seed namespace에서 실행한다. 기존 MM-main의 SCM1,
$n_R=100$ shared·shifted 결과는 S01·S02의 replication으로 재사용하지 않는다.

### 5.3 Cell과 route의 구분

Covariate-shift cell에서는 같은 원자료에 다음 세 ratio route를 적용한다.

1. `exact`: analytic $r_0$
2. `classifier`: nuisance 역할에서 적합하고 이후 freeze한 density-ratio estimator
3. `wrong_unit`: $r(x)=1$로 두는 negative control

세 route는 세 개의 다른 데이터 cell이 아니다. 같은 `cell_id`, `rep_id`,
`bundle_id` 아래 생성된 데이터에 적용되는 세 분석 경로다. 따라서 6개
covariate-shift cell이 18개 cell로 늘어나지 않는다.

Shared와 outcome-mismatch cell의 primary ratio route는 `exact_unit`이다.
Non-transport comparator를 추가하더라도 replication 수로 세지 않는다.

---

## 6. 동일한 ATE-CATE role bundle

### 6.1 표본 역할

한 replication은 다음 네 역할을 독립적으로 생성한다.

| role | 용도 | RCT 크기 | OBS 크기 |
|---|---|---:|---:|
| nuisance | outcome regression, basis, density ratio 적합 | $n_R$ | 3,000 |
| tuning | ATE coefficient와 CATE candidate 적합 | $n_R$ | 1,000 |
| selection | CATE candidate 선택 | $n_R$ | 1,000 |
| evaluation/reporting | 최종 ATE 평가와 CATE reporting score | $n_R$ | 1,000 |

OBS 크기 벡터는 모든 synthetic cell에서

$$
(N_O^{\mathrm{nuis}},N_O^{\mathrm{tune}},N_O^{\mathrm{select}},
N_O^{\mathrm{eval}})=(3000,1000,1000,1000)
$$

으로 고정한다.

### 6.2 동일 bundle의 의미

같은 `(protocol_id, phase, cell_id, rep_id)`에서 ATE와 CATE는 다음을 공유한다.

- source-law parameter
- role별 row draw
- potential outcome noise
- treatment assignment
- nuisance fits와 density-ratio fits, estimator 정의가 허용하는 범위

각 결과 행에는 같은 `bundle_id`와 role fingerprint를 기록한다. ATE와 CATE의
비교는 이렇게 paired된다.

### 6.3 독립성

공유는 역할 사이의 재사용을 뜻하지 않는다. seed key는 최소한 다음을 포함한다.

```text
protocol_id / phase / cell_id / rep_id / source / role / component / substep
```

Nuisance, tuning, selection, evaluation/reporting stream은 서로 달라야 한다.
Evaluation outcome만 바꿨을 때 nuisance fit, candidate grid와 selected candidate가
변하지 않아야 한다.

### 6.4 Truth integration

Synthetic CATE risk는 각 SCM의 RCT target law에서 생성한 50,000-point 독립 truth
grid로 적분한다. SCM마다 하나의 grid를 main performance seed와 다른 namespace에서
생성해 freeze하고 selection에는 제공하지 않는다. Grid 크기, 생성법 또는 numerical
integration tolerance를 바꾸려면 smoke 전에 protocol version을 바꿔야 한다.
Smoke 이후에는 바꾸지 않는다.

---

## 7. 세 regime의 생성 계약

### 7.1 Shared nominal

Shared cell에서는

$$
P_R^X=P_O^X,\qquad r_0(x)=1
$$

이다. RCT와 OBS의 outcome mechanism은 해당 SCM의 nominal 정의를 따른다.
`exact_unit` route가 primary다.

### 7.2 Covariate-shift nominal

SCM $k$의 frozen nominal trial law와 observational law를 다음처럼 둔다.

$$
X_R\sim\mathcal N(\mu_k,\Sigma_k),\qquad
X_O\sim\mathcal N(\mu_k+d_k,\Sigma_k).
$$

기존 SCM별 shift direction을 $v_k$라 하고

$$
d_k=
\frac{v_k}{\sqrt{v_k^\top\Sigma_k^{-1}v_k}}
$$

로 고정한다. 따라서 shift의 squared Mahalanobis distance는

$$
d_k^\top\Sigma_k^{-1}d_k=1
$$

이다. Exact target-to-source ratio는

$$
r_{0,k}(x)
=\exp\left[
-d_k^\top\Sigma_k^{-1}(x-\mu_k)
+\frac12d_k^\top\Sigma_k^{-1}d_k
\right].
$$

이 shift는 $P_R^X\ne P_O^X$와 $P_R^X\ll P_O^X$를 만족한다. $v_k$, $d_k$와
analytic ratio formula는 smoke 전에 fingerprint에 고정하고 smoke 또는 pilot 뒤에
바꾸지 않는다. Outcome law와 true $\tau(x)$는 shared nominal cell과 동일하게
유지한다.

같은 `(SCM, n_R, rep_id)`의 shared·shifted pair에는 가능한 모든 base randomness에
common random numbers를 사용한다. 바뀌는 것은 protocol에 선언된 covariate-law
transformation뿐이다. Pair fingerprint가 이를 증명해야 한다.

세 ratio route는 동일한 shifted data를 받는다. `wrong_unit` 결과가 나쁘지 않아도
engineering failure가 아니며, 나쁘다고 해서 exact 또는 classifier route의
validity가 자동으로 증명되는 것도 아니다.

### 7.3 Pure outcome-domain mismatch

이 cell은 covariate shift를 제거한다.

$$
P_R^X=P_O^X,\qquad r_0(x)=1.
$$

SCM $k\in\{1,2,3\}$에 대해 trial covariate law의 첫 좌표 평균과 분산을
$\mu_{k1}$과 $\Sigma_{k,11}$로 쓰고

$$
h_k(X)=\frac{X_1-\mu_{k1}}{\sqrt{\Sigma_{k,11}}}
$$

를 정의한다. 따라서 population law에서

$$
\mathbb E_R\{h_k(X)\}=0,\qquad
\operatorname{Var}_R\{h_k(X)\}=1.
$$

OBS의 두 potential outcome에 같은 baseline perturbation을 더한다.

$$
Y_O(a)=Y_O^{\mathrm{base}}(a)+\delta_k h_k(X),
\qquad a\in\{0,1\}.
$$

강도는

$$
\delta_k
=\texttt{SCMParameters.outcome\_noise\_sd}
=\sigma_{\varepsilon,k}
$$

로 고정한다. 즉 $\delta_k$는 SCM $k$에 이미 정의된 outcome noise standard
deviation $\sigma_{\varepsilon,k}$ 자체다. Arm mixture, 표본 standard deviation 또는
별도 calibration estimate로 다시 정의하지 않는다. Smoke와 pilot 결과를 보고
바꾸지 않는다.

Perturbation은 두 arm에 동일하므로

$$
\begin{aligned}
\tau_O(X)
&=Y_O(1)-Y_O(0)\\
&=Y_O^{\mathrm{base}}(1)-Y_O^{\mathrm{base}}(0)
\end{aligned}
$$

가 성립해야 한다. Raw potential-outcome fixture에서 machine precision으로 검사한다.

$h_k$는 $X_1$의 affine transformation이므로 frozen OBS linear feature span에
정확히 들어가야 한다. Nuisance sample을 보기 전에 feature schema를 정하고,
basis는 nuisance 단계 이후 freeze하여 tuning, selection, evaluation에서
재추정하지 않는다. Shared와 outcome-mismatch pair도 같은 `(SCM,n_R,rep_id)`에서
공변량, treatment와 base outcome noise에 common random numbers를 사용한다.
두 cell 사이에서 바뀌는 것은 위의 OBS perturbation뿐이다. 이 설계는
covariate imbalance가 아니라 outcome-domain mismatch에 대한 반응을 분리한다.

---

## 8. Estimator와 비교 경로

### 8.1 CATE learner와 candidate

CATE는 DRF와 RF를 각각 평가한다. 두 learner 모두 fixed cubic spline
$K_{\mathrm{spline}}=3$을 사용한다. Neural candidate는 이 PRD의 primary design에
포함하지 않는다.

CATE candidate grid는

$$
(\lambda,\omega)\in
\{0,0.25,0.5,0.75,1\}^2
$$

로 고정한다. 각 restriction에 맞는 grid subset 안에서 selection sample로
candidate를 고른다.

| label | coefficient restriction | 의미 |
|---|---|---|
| `rct_only` | $(\lambda,\omega)=(0,0)$ | RCT-only 기준 |
| `lambda_only` | $\omega=0$ | outcome-regression blend만 허용 |
| `omega_only` | $\lambda=0$ | transported control variate만 허용 |
| `joint` | $\lambda,\omega$ 모두 선택 | 두 channel을 함께 허용 |
| `cate_grid_oracle` | true CATE risk로 grid에서 선택 | synthetic CATE diagnostic only |

`cate_grid_oracle`은 실행 가능한 estimator가 아니며 selection 결과에 섞지 않는다.
`selected`와 `cate_grid_oracle`은 모든 raw·summary·plot에서 다른 `method_label`을
쓴다. CATE selection regret도 이 grid oracle을 기준으로만 정의한다.

### 8.2 ATE coefficient와 oracle

ATE는 CATE의 discrete grid selection을 사용하지 않는다. Tuning 역할에서 ATE의
plug-in variance criterion을 구성한 뒤 각 restriction 아래

$$
(\lambda,\omega)\in[0,1]^2
$$

로 projection한 minimizer를 사용한다.

- `rct_only`: $(0,0)$
- `lambda_only`: $\omega=0$에서 projected plug-in minimization
- `omega_only`: $\lambda=0$에서 projected plug-in minimization
- `joint`: $[0,1]^2$에서 projected plug-in minimization

Population 또는 별도 large-Monte-Carlo moments에서 계산한 ATE oracle은
`ate_population_oracle`로 따로 보고한다. 이것은 구현 가능한 plug-in estimator와
다른 diagnostic이다. ATE에는 `grid_oracle` 또는 CATE selection regret를 적용하지
않는다.

### 8.3 Same-$\lambda$ ablation

Joint-selected CATE candidate의 $\lambda$를 고정하고 같은 grid에서 $\omega=0$인
기존 candidate와 비교하는 ablation은 secondary diagnostic으로만 허용한다.

다음 조건을 모두 만족할 때만 계산한다.

1. 비교점이 원래 사전 고정 candidate grid에 들어 있다.
2. 새 적합이나 새로운 outcome inspection이 필요하지 않다.
3. primary `lambda_only`와 혼동되지 않는 label을 쓴다.
4. confirmatory success criterion에 넣지 않는다.

조건을 만족하지 않으면 생략한다.

ATE에서는 tuning 역할에서 얻어 freeze한 joint coefficient의 첫 좌표를
$\widehat\lambda_{\mathrm{joint}}$라 한다. 같은 reporting/evaluation bundle에서
재적합이나 재선택 없이

$$
\widehat\theta_r(\widehat\lambda_{\mathrm{joint}},0)
$$

을 평가한다. 이 estimate의 squared error와 같은 reporting/evaluation 값으로
계산한 대응 variance estimate를 모두
`ate_joint_lambda_omega0_secondary`라는 별도 label 아래 기록한다. Primary
`lambda_only` 또는 `joint`를 덮어쓰지 않고 confirmatory success criterion에도
넣지 않는다.

### 8.4 Learner scope

DRF와 RF는 각각 원고에 정의된 loss와 score를 사용한다. Primary CATE 구현은 위의
fixed spline만 사용한다. 이 PRD는 neural empirical extension이나 새 theorem을
추가하지 않는다.

### 8.5 Clipping과 numerical failure

다음을 지킨다.

- outcome prediction과 density ratio를 임의로 자르지 않는다.
- spline basis boundary 처리는 outcome 또는 ratio clipping과 별도로 기록한다.
- nonfinite, singular solve, optimizer failure를 임의의 유한값으로 바꾸지 않는다.
- fallback estimator를 사용하려면 사전에 method label과 조건을 protocol에 넣는다.
- 실패한 replication을 평균에서 조용히 제외하지 않는다.

---

## 9. 실행 단계와 계산 장부

### 9.1 단계 원칙

Smoke와 pilot은 engineering·runtime·tail 진단용 exploratory 자료다. Main은 별도
confirmatory namespace에서 시작한다. 세 phase 사이에는 performance row를 합치지
않는다.

### 9.2 Smoke

| 항목 | 계약 |
|---|---|
| synthetic | 18 cells × 1 rep = 18 cell-replications |
| STAR | 1 conditional bundle |
| 목적 | schema, seed, role separation, data identity, finite path, failure serialization |
| scientific inference | 금지 |

Smoke PASS는 모든 hard assertion이 통과하고 모든 scheduled unit이 성공 또는
명시적 failure로 기록됐다는 뜻이다. Method가 RCT보다 좋아야 한다는 조건은 없다.

### 9.3 Pilot

| 항목 | 계약 |
|---|---|
| synthetic | 18 cells × 20 reps = 360 cell-replications |
| STAR | 20 conditional bundles |
| 목적 | runtime, memory, failure rate, tail, aggregation과 plot pipeline 점검 |
| main pooling | 금지 |

Smoke 전에 candidate grid, $h_k$, $\delta_k$, shift direction, ratio route, metric,
uncertainty unit과 truth grid를 freeze한다. Primary 95% interval은 whole paired
bundle을 재추출하는 percentile bootstrap으로 계산하며, resample 수는 $B=2000$이다.
Bootstrap RNG는 protocol·analysis·estimand·cell 또는 predeclared aggregate key의
SHA256에서 유도한 seed를 사용한다. Pilot outcome을 보고 이것을 바꾸려면 protocol을
새 버전으로 만들고 기존 pilot은 개발 자료로만 남긴다.

### 9.4 Main

| 항목 | 계약 |
|---|---|
| synthetic | 18 cells × 200 reps = 3,600 cell-replications |
| STAR | 총 200 conditional bundles |
| seed | smoke·pilot과 disjoint confirmatory namespace |
| extension | 없음 |

Main의 $K=200$은 고정이다. 결과가 유리하거나 불리하거나 불확실하다는 이유로
연장하지 않는다. 추가 연구가 필요하면 별도 protocol과 별도 claim family로 한다.

### 9.5 무엇을 replication으로 세는가

`cell_id`와 `rep_id`로 식별되는 source-data bundle 한 개가 한 replication이다.
같은 bundle에서 나온 learner, method, ratio route, metric 행은 replication이 아니다.

Main synthetic 장부는 정확히 다음과 같다.

$$
18\text{ cells}\times200\text{ replications}=3{,}600
\text{ cell-replications}.
$$

운영상 main source-data bundle은

$$
3{,}600\text{ synthetic}+200\text{ conditional STAR}=3{,}800
$$

개다. 이 합은 실행·storage·failure accounting에만 사용한다. Synthetic과 STAR는
서로 다른 반복 표본과 estimand를 가지므로 effect estimate, variance, interval 또는
유효 표본 크기에서 절대 pooling하지 않는다.

---

## 10. Synthetic 평가량

### 10.1 ATE 지표

Replication $s$의 ATE estimate를 $\widehat\theta_s$, 참값을 $\theta$라 하면 error는

$$
e_s=\widehat\theta_s-\theta
$$

이다. Cell-level 지표는 다음과 같다.

$$
\operatorname{Bias}=\frac1K\sum_{s=1}^K e_s,
$$

$$
\operatorname{Variance}=\frac1{K-1}\sum_{s=1}^K
(\widehat\theta_s-\overline{\widehat\theta})^2,
$$

$$
\operatorname{MSE}=\frac1K\sum_{s=1}^K e_s^2,
\qquad
\operatorname{RMSE}=\sqrt{\operatorname{MSE}}.
$$

Interval estimator가 protocol대로 구현된 method에는 empirical coverage와 평균
interval length를 보고한다. 구현하지 않은 method에는 `N/A`를 기록하고 임의의
interval을 사후 추가하지 않는다.

Primary paired contrast는 같은 bundle의 joint와 RCT-only squared error 차이다.

$$
\Delta^{\mathrm{ATE}}_s
=
(\widehat\theta_{s,\mathrm{joint}}-\theta)^2
-(\widehat\theta_{s,\mathrm{RCT}}-\theta)^2.
$$

음수는 joint의 squared error가 더 작다는 뜻이다.

### 10.2 CATE 지표

같은 frozen target grid에서 replication $s$의 true integrated squared risk는

$$
R_s(\widehat\tau)
=
\mathbb E_R
\left[
\{\widehat\tau_s(X)-\tau(X)\}^2
\right]
$$

이다. Numerical grid average를 population expectation의 근사로 사용하고 grid
fingerprint와 integration error diagnostic을 기록한다.

Primary paired contrast는

$$
\Delta^{\mathrm{CATE}}_s
=
R_s(\widehat\tau_{\mathrm{joint}})
-R_s(\widehat\tau_{\mathrm{RCT}})
$$

이다. 음수는 joint risk가 더 작다는 뜻이다.

Risk ratio도 보고한다.

$$
Q_s^{\mathrm{CATE}}
=
\frac{R_s(\widehat\tau_{\mathrm{joint}})}
{R_s(\widehat\tau_{\mathrm{RCT}})}.
$$

분모가 너무 작거나 nonfinite이면 ratio를 강제로 만들지 않고 status를 기록한다.
Primary inference는 additive risk difference를 사용한다.

### 10.3 CATE variance와 squared bias

고정 target point $x$에서

$$
\overline{\widehat\tau}(x)=\frac1K\sum_{s=1}^K\widehat\tau_s(x)
$$

로 두고 finite-$K$ prediction variance를

$$
\widehat V(x)=\frac1K\sum_{s=1}^K
\left\{
\widehat\tau_s(x)-\overline{\widehat\tau}(x)
\right\}^2
$$

로 계산한다. Squared bias는

$$
\widehat B^2(x)=
\left\{
\overline{\widehat\tau}(x)-\tau(x)
\right\}^2
$$

이다. Target grid에서 각각 평균하고, 같은 finite-$K$ convention을 모든 method에
적용한다. Risk, variance, squared bias의 numerical identity residual도 저장한다.

### 10.4 CATE selection regret

Selection sample로 고른 candidate를 $\widehat c_s$, 동일한 사전 고정 grid에서
true risk가 가장 작은 candidate를 $c_s^*$라 하면

$$
\operatorname{Regret}_s
=R_s(\widehat c_s)-R_s(c_s^*)
$$

이다. $c_s^*$는 CATE 평가용 oracle label이며 selection 과정에 들어가지 않는다.
이 regret는 ATE나 STAR에 적용하지 않는다.

### 10.5 Ratio와 tail 진단

Applicable route에는 다음을 저장한다.

- normalization error
- balance residual
- effective sample size와 source-size 대비 비율
- maximum weight
- nonfinite weight count
- clipping fraction, protocol상 clipping이 없으면 정확히 0
- optimizer success, protocol status, raw message, iterations와 objective
- CATE risk와 ATE squared error의 median, 90th, 95th, 99th percentile와 maximum
- paired RCT보다 나쁜 replication 비율
- unit별 wall-clock runtime과 peak memory가 가능하면 그 값

---

## 11. 불확실성 단위와 판정

### 11.1 독립 단위

Primary interval은 95% percentile paired bootstrap이다. Resampling unit은 whole
replication bundle이고 resample 수는 정확히 $B=2000$이다. 한 bundle 안의 ATE,
CATE, learner, method, ratio route를 따로 resample하지 않는다. 각 bootstrap seed는
protocol fingerprint와 analysis key의 hash로 결정되며 artifact에 기록한다.

Cell-level 결과는 해당 cell의 200 paired bundle을 재추출한다. SCM 전체 또는
regime 전체처럼 protocol에서 미리 지정한 aggregate만 SCM family 안에서
family-stratified paired bootstrap을 사용한다. 지정되지 않은 사후 aggregate에는
이를 적용하지 않는다.

### 11.2 Engineering PASS

다음이 모두 충족되면 결과 방향과 무관하게 engineering PASS다.

1. 18개 cell과 cell당 200개의 main seed가 정확히 존재한다.
2. 모든 seed는 result 또는 explicit failure를 가진다.
3. Role disjointness와 ATE-CATE bundle identity가 검증된다.
4. Shifted 세 route가 같은 raw-data fingerprint를 가진다.
5. Outcome mismatch의 표준화와 treatment-effect invariance가 검증된다.
6. 모든 fit·selection·evaluation은 정해진 역할만 사용한다.
7. Hidden clipping, redraw, dropped failure가 없다.
8. Raw에서 summary, 표, 그림을 재생성한 값이 일치한다.

Joint가 RCT-only보다 나쁜 것은 engineering failure가 아니다.

### 11.3 과학적 방향 label

Loss difference의 사전 고정 uncertainty interval을 $[L,U]$라 한다. 방향 label은
다음처럼 기계적으로 붙인다.

| 조건 | label | 허용 해석 |
|---|---|---|
| $U<0$ | `COMPATIBLE_WITH_IMPROVEMENT` | 평균 loss 감소와 양립 |
| $L\le0\le U$ | `INCONCLUSIVE_DIRECTION` | 방향 불확실 |
| $L>0$ | `COMPATIBLE_WITH_HARM` | 평균 loss 증가와 양립 |

이 label은 descriptive interval label이며 multiplicity-adjusted hypothesis decision이
아니다. Confirmatory multiplicity family와 실질적 threshold를 사용할 경우 smoke와
pilot outcome을 보기 전에 별도 decision record에 고정한다. 그런 record가 없으면
방향 label만 보고하고 familywise 또는 confirmatory significance 주장을 하지 않는다.

### 11.4 Mechanism interpretation

- shared 대 covariate shift: transport sensitivity
- exact 대 classifier: density-ratio estimation cost
- exact 대 wrong-unit: transport misspecification negative control
- shared 대 outcome mismatch: OBS baseline compatibility sensitivity
- `lambda_only`, `omega_only`, `joint`: 두 channel의 구분 가능한 기여

한 비교의 결과로 다른 비교의 원인을 확정하지 않는다. 예를 들어 classifier와
exact가 비슷하다는 사실은 일반적으로 ratio estimation이 공짜라는 뜻이 아니다.

---

## 12. STAR real-outcome stress test

### 12.1 입력과 frozen source law

STAR primary outcome은 `mathk` 하나다. Treatment는 small class 대 pooled
regular/regular+aide다.

입력은 `--star-csv` 또는 `DATAFUSIONPPI_STAR_CSV`로 지정된 기존 로컬 파일을
해결한다. Scientific run 중 다운로드하거나 project로 복사하지 않는다. 현재
[cate_star.py](../code/confirm/cate_star.py#L21-L55)는 source URL, SHA256, row count와
필수 schema 검사를 제공한다.

Performance bundle 전에 다음을 한 번 freeze한다.

- eligible exact-X cohort
- ordered row ID hash
- OOF fold assignment hash
- exact-X probability와 propensity table hash
- source weight와 outcome-standardization metadata
- learner feature tuple과 missingness indicators

현재 exact-X construction과 fingerprint fields는
[cate_star.py](../code/confirm/cate_star.py#L119-L217)를 근거로 한다.

### 12.2 Source law의 의미

Frozen law는 source 수준에서

$$
P_O^X=P_R^X,\qquad r_0(x)=1
$$

을 강제한다. Exact-cell trial propensity는 eligible cell에서

$$
0.15\le e_x\le0.85
$$

를 만족해야 한다. Realized finite sample의 covariate share가 정확히 같을 필요는
없다.

이 구성은 source marginal equality를 보장한다. STAR treatment가 모든 potential
outcome에 대해 교환 가능하다는 것, 또는 causal truth를 안다는 것을 보장하지
않는다.

### 12.3 Conditional resplit bundle

STAR main은 200개의 독립 conditional resplit/resample bundle을 사용한다.

- cohort와 source-law freeze는 200회에 공통이다.
- 각 bundle의 nuisance, tuning, selection, evaluation/reporting draw는 독립이다.
- ATE와 CATE는 같은 bundle ID와 role draws를 공유한다.
- sampling은 frozen empirical support에 조건부인 with-replacement draw다.
- $K=200$에서 종료하며 결과에 따라 연장하지 않는다.

따라서 200회는 새로운 학교 모집단을 200번 독립 추출한 반복이 아니다. 고정된
eligible cohort와 frozen empirical source law 아래의 conditional resampling
variability를 나타낸다.

### 12.4 STAR에서 허용하는 지표

다음만 primary 또는 supporting metric으로 허용한다.

- 원고의 ATE estimator가 해당 STAR construction에 적용되는 각 method의 descriptive
  ATE estimate
- 각 bundle에서 원고 variance formula가 적용되는 method의 plug-in variance,
  standard error와 95% interval
- 200개 estimate의 across-resplit empirical variance와 stability summary
- 동일 bundle에서 계산한 method 간 descriptive ATE difference
- Independent evaluation/reporting sample의 CATE risk-difference score
- $\lambda,\omega$ 선택 분포
- exact-X propensity, ratio, ESS, support와 role diagnostics
- failure rate와 runtime

ATE estimate, variance, standard error, interval과 across-resplit stability에는 모두
`conditional on the frozen eligible STAR cohort and source law` label을 붙인다.
원고 estimator의 가정이나 variance formula가 적용되지 않는 method에는 SE·CI를
만들지 않고 `N/A`와 이유를 기록한다. Across-resplit variation은 모집단 sampling
variance 또는 per-bundle plug-in variance와 같은 양으로 부르지 않는다. STAR의
CATE 성능량은 reporting-score difference 하나뿐이다. Synthetic true-risk, bias²,
variance decomposition 또는 oracle regret를 대신 붙이지 않는다.

### 12.5 STAR에서 금지하는 지표와 표현

다음은 보고하지 않는다.

- true CATE risk
- true CATE squared bias
- true treatment-effect coverage
- synthetic truth에 대한 RMSE
- 모집단 causal ground truth를 회복했다는 주장
- source-law equality만으로 causal identification이 성립한다는 주장
- conditional resplit interval을 모집단 sampling interval로 해석하는 문장

STAR held-out score는 synthetic true risk와 다른 quantity다. 같은 축이나 같은
평균으로 합치지 않는다.

---

## 13. MM-main과 기존 ATE 결과의 경계

### 13.1 MM-main은 external consistency only

[MM-main RECORD](confirm/MM-main-66e93079cbd7/RECORD.md)는 SCM1, 역할별
$n_R=100$의 shared와 shifted 결과를 제공한다. 새 S01·S02와 비교할 때 다음만
허용한다.

- estimator와 DGP mapping이 일치하는지 확인
- direction과 effect magnitude가 명백하게 모순되는지 기술
- 코드 또는 protocol 변화가 차이를 설명할 수 있는지 진단

다음은 금지한다.

- 기존 400 cell record를 새 3,600 record에 합치기
- 기존 bootstrap draw를 새 interval에 넣기
- 기존 결과를 새 seed prefix로 간주하기
- 새 결과가 불리할 때 기존 결과를 선택적으로 보강 자료로 합산하기

### 13.2 기존 ATE는 별도 evidence

기존 ATE 결과는 구현과 예상 규모를 확인하는 배경자료다. 새 18-cell 분석은
동일 ATE-CATE bundle이라는 새 설계 목적이 있으므로 ATE도 main에서 새로 계산한다.
기존 ATE estimate를 새 CATE bundle과 인위적으로 짝짓지 않는다.

---

## 14. 필수 산출물 계약

최종 hard output은 그림 3개와 표 2개다. 추가 부록은 가능하지만 이 다섯 개를
대체하지 않는다.

### Figure 1. Synthetic primary performance

목적은 18-cell rectangle에서 ATE와 CATE의 primary paired contrast를 보여주는
것이다.

- Panel A1: ATE joint-minus-RCT paired squared-error difference
- Panel A2: joint와 `rct_only`의 empirical ATE variance
- Panel B: CATE joint-minus-RCT paired true-risk difference
- 행 또는 facet: SCM1, SCM2, SCM3
- 열 또는 shape: $n_R=100,400$
- color: shared, covariate shift, outcome mismatch
- Panel A1·B point: bundle-level loss contrast의 평균
- Panel A2 point: $K_{\mathrm{rep}}=200$개 ATE estimate에 걸친 empirical variance
- interval: 사전 고정 paired uncertainty interval
- vertical reference: 0

Panel A1과 B는 bundle마다 먼저 loss contrast를 만든 뒤 그 200개 값의 평균을
표시한다. Panel A2는 per-bundle variance estimate의 평균이 아니다. 각 bootstrap
draw에서 whole bundles를 재추출하고, 그 draw 안의 200개 ATE estimate로 joint와
`rct_only`의 empirical variance를 각각 다시 계산한다. A2에는 matched-RCT variance
ratio를 보조 숫자로 주석 처리할 수 있지만 primary point는 두 absolute variance다.

ATE squared-error difference, ATE empirical variance, CATE risk difference의 y축은
각각 분리한다. 서로 다른 metric을 표준화하지 않은 채 하나의 numerical effect로
합치지 않는다.

### Figure 2. Channel과 failure mechanism

목적은 $\lambda$와 $\omega$ channel 및 transport failure를 분해하는 것이다.

- `rct_only`, `lambda_only`, `omega_only`, `joint`의 relative performance
- CATE risk의 finite-$K$ variance와 squared-bias decomposition
- covariate-shift cell의 exact, classifier, wrong-unit route
- outcome-mismatch cell의 coefficient 선택률
- tail 또는 failure marker가 있는 경우 평균만으로 숨기지 않음

Same-$\lambda$, $\omega=0$ ablation을 표시하면 secondary임을 figure 안에 적는다.

### Figure 3. STAR conditional resplit result

- Panel A: descriptive ATE의 method 차이와 conditional interval
- Panel B: held-out CATE reporting-score difference
- Panel C: $\lambda,\omega$, ESS, failure 또는 support diagnostic

Figure title과 caption에 `conditional on the frozen eligible STAR cohort and source
law`를 명시한다. True CATE risk라는 label을 쓰지 않는다.

### Table 1. Synthetic confirmatory results

정확히 18개 cell을 기본 행 단위로 하며 다음을 포함한다.

- SCM, $n_R$, regime, ratio route, learner, method
- 독립 replication 수와 failure 수
- ATE bias, joint와 `rct_only`를 포함한 method별 absolute variance, matched-RCT
  variance ratio, MSE, RMSE, coverage와 interval length
- CATE risk, prediction variance, squared bias, regret
- paired contrast와 uncertainty interval
- ratio와 coefficient 핵심 진단

Route와 method 때문에 행이 늘어나도 `independent_replications=200`을 명시한다.

### Table 2. STAR, robustness와 provenance

다음을 포함한다.

- STAR descriptive ATE와 held-out score
- conditional bundle 수와 failure 수
- selected coefficient 분포
- exact-X propensity, support, ESS와 source-law diagnostic
- code, protocol, dataset, cohort, seed와 artifact fingerprint
- MM-main external consistency comparison, 새 추정치와 분리된 block

MM-main을 pooled sample size에 포함하지 않는다.

---

## 15. Claim-to-evidence map

Claim을 만들 때 다음 literal cell mapping만 사용한다.

- Covariate-shift evidence: S02, S05, S08, S11, S14, S17
- Outcome-mismatch evidence: S03, S06, S09, S12, S15, S18
- Shared 대 outcome-mismatch paired cells:
  S01↔S03, S04↔S06, S07↔S09, S10↔S12, S13↔S15, S16↔S18
- $n_R$ comparison은 같은 SCM과 regime 안에서만 한다:
  S01↔S04, S02↔S05, S03↔S06;
  S07↔S10, S08↔S11, S09↔S12;
  S13↔S16, S14↔S17, S15↔S18
- STAR evidence는 정확히 200 conditional bundles뿐이며 synthetic cell mapping에
  넣지 않는다.

| 주장 후보 | 필요한 근거 | 허용 범위 | 금지되는 확대 |
|---|---|---|---|
| 특정 SCM에만 의존하지 않는다 | 세 SCM에서 같은 방향과 적절한 uncertainty | 이 세 SCM family | 모든 DGP에 일반화 |
| 작은 trial에만 국한되지 않는다 | $n_R=100,400$ 비교 | 두 역할별 표본 크기 | 임의의 $n_R$ |
| exact transport에서 target을 유지한다 | exact route의 ATE bias·CATE risk와 구현 identity | 실행한 support와 DGP | estimated ratio도 exact라는 주장 |
| classifier 비용이 작다 | paired exact 대 classifier | 사용한 classifier와 shift | ratio estimation이 일반적으로 공짜 |
| wrong transport는 위험할 수 있다 | wrong-unit negative control의 tail과 평균 | 지정한 shift | 항상 harmful |
| 두 channel은 구분되는 기회를 가진다 | shared·mismatch 및 channel ablation | 사전 고정 grid | structural identification의 증명 |
| outcome mismatch에 강건하다 | 6개 mismatch cell의 ATE·CATE 결과 | $\delta_k=\sigma_{\varepsilon,k}$와 지정한 $h_k$ | 임의 outcome drift |
| ATE와 CATE가 함께 개선된다 | 동일 bundle의 분리된 두 paired contrast | 두 estimand에 대한 병렬 evidence | 두 estimand가 동일 |
| STAR에서 유용하다 | descriptive ATE와 held-out score | frozen cohort 조건부 | causal truth 회복 |
| 원고 theorem과 양립한다 | theorem eligibility audit | 해당 exact assumptions와 learner | empirical 결과를 proof로 표현 |

원고의 현재 실험 서술은 CATE가 계획 단계임을 밝힌다. 이 PRD의 실행이 끝나도
원고를 자동 변경하지 않는다. Claim map과 최종 evidence를 별도 검토한 뒤 원고
수정을 결정한다.

---

## 16. Artifact와 provenance 계약

### 16.1 권장 stage tree

아래는 planned artifact contract다. 실제 경로 이름에는 protocol fingerprint를
포함한다.

```text
materials/confirm/<stage>-<fingerprint>/
  FINGERPRINT
  RECORD.md
  SHA256SUMS
  code_manifest.json
  protocol.json
  schedule.json
  seed_ledger.json
  environment.json
  input_manifest.json
  role_fingerprints.json
  checks.json
  failures.json
  runtime.json
  raw/
    ate_results.csv
    cate_results.csv
    coefficient_results.csv
    ratio_diagnostics.csv
    sufficient_statistics.json
  analysis/
    summary.json
    table1_synthetic.csv
    table2_star_provenance.csv
    figure1_synthetic_primary.png
    figure2_mechanisms.png
    figure3_star.png
    analysis_provenance.json
```

STAR와 synthetic을 별도 stage로 두더라도 최종 manifest는 양쪽 fingerprint를
참조해야 한다.

### 16.2 Fingerprint

Fingerprint에는 최소한 다음을 넣는다.

- protocol 문서와 machine-readable config hash
- 실행 code와 requirements hash
- DGP와 candidate grid
- 모든 sample-role size
- $h_k$, $\mu_{k1}$, $\Sigma_{k,11}$, $\delta_k$, $\sigma_{\varepsilon,k}$ 정의
- truth-grid config와 hash
- STAR input·cohort·fold·probability·propensity hash
- seed-key schema와 schedule hash
- 분석 metric과 uncertainty specification

### 16.3 Immutable raw와 no overwrite

- Existing promoted stage는 덮어쓰지 않는다.
- 같은 fingerprint의 완전한 stage가 있으면 read-only validate한다.
- 불완전 stage를 재개할 때도 기존 success row를 수정하지 않는다.
- 새 run은 새 scratch directory에 쓴다.
- 모든 required output을 검증한 뒤에만 atomic promotion한다.
- Promotion 실패 시 새 run이 소유한 scratch만 정리한다.
- 기존 scratch나 historical stage를 scan-delete하지 않는다.

### 16.4 Failure와 resume

모든 scheduled seed에는 다음 중 하나가 있어야 한다.

1. 완전한 result record
2. stage, cell, rep, component, exception type과 message를 가진 failure record

Failure seed를 새 seed로 바꾸지 않는다. Resume은 ledger에서 완전하고 hash가 맞는
unit만 건너뛴다. 부분 output은 성공으로 간주하지 않는다.

### 16.5 결정적 산출물

Timestamp와 runtime을 제외한 scientific artifact는 같은 code, protocol, input,
seed에서 byte-identical해야 한다. 작은 fixed-seed fixture로 이를 검증한다. 전체
main을 두 번 실행하는 것을 기본 determinism test로 요구하지 않는다.

---

## 17. Acceptance-test matrix

| ID | 검사 | 대상 | PASS 기준 | 실패 시 조치 |
|---|---|---|---|---|
| G01 | cell completeness | synthetic schedule | 정확히 18 unique cells | 실행 금지 |
| G02 | Cartesian product | synthetic schedule | SCM 3 × $n_R$ 2 × regime 3 | 실행 금지 |
| G03 | role size·use | all synthetic cells | RCT 각 $n_R$, OBS 3000/1000/1000/1000이고 nuisance가 ratio를 적합 | 실행 금지 |
| G04 | role disjointness | each bundle | row와 seed overlap 0 | stage failure |
| G05 | ATE-CATE identity | each bundle | source·role fingerprint 동일 | stage failure |
| G06 | route pairing | shifted cells | 세 route raw-data hash 동일 | stage failure |
| G07 | exact unit ratio | shared·mismatch | $r_0=1$ identity | stage failure |
| G08 | shift ratio identity | shifted exact | Mahalanobis distance 1과 analytic density-ratio fixture 통과 | stage failure |
| G09 | classifier diagnostics | shifted classifier | convergence·objective·weight diagnostics 존재 | failure 기록 |
| G10 | no hidden clipping | every route | clip fraction 0, protocol과 code 일치 | stage failure |
| G11 | $h_k$ centering | mismatch | analytic population mean 0과 fixture 일치 | 실행 금지 |
| G12 | $h_k$ scaling | mismatch | analytic population variance 1과 fixture 일치 | 실행 금지 |
| G13 | feature inclusion | mismatch | frozen feature schema에 $h_k$ 존재 | stage failure |
| G14 | effect invariance·CRN | mismatch | 양 arm perturbation 차이 0이고 shared pair의 base fingerprint 일치 | stage failure |
| G15 | selection honesty | all | reporting perturbation이 selected candidate를 변경하지 않음 | stage failure |
| G16 | oracle separation | synthetic | CATE selected·grid oracle와 ATE plug-in·population oracle 분리 | stage failure |
| G17 | seed accounting | every phase | scheduled = success + failure | promotion 금지 |
| G18 | no redraw | every phase | failure seed가 ledger에 그대로 존재 | promotion 금지 |
| G19 | main count | synthetic | cell당 exactly 200, 총 3,600 | promotion 금지 |
| G20 | no extension | main | rep ID가 fixed schedule 밖에 없음 | promotion 금지 |
| G21 | phase separation | all | smoke·pilot·main seed namespace disjoint | promotion 금지 |
| G22 | ATE identity | fixed arrays | estimator와 direct formula 일치 | 실행 금지 |
| G23 | ATE variance | fixed arrays | production과 independent formula 일치 | 실행 금지 |
| G24 | CATE risk | fixed grid | raw prediction과 summary risk 일치 | promotion 금지 |
| G25 | decomposition | CATE summary | risk와 variance+bias² identity residual 허용치 이내 | promotion 금지 |
| G26 | paired bootstrap | analysis | 95% percentile, $B=2000$, hash seed, whole bundle resampling | promotion 금지 |
| G27 | failure tails | analysis | quantiles가 failure omission 없이 계산됨 | promotion 금지 |
| G28 | MM exclusion | primary analysis | MM-main row count 0 | promotion 금지 |
| G29 | STAR input | STAR preflight | file, SHA256, row count, schema 일치 | STAR 실행 금지 |
| G30 | STAR freeze | STAR | support·fold·probability·propensity hash 존재 | STAR 실행 금지 |
| G31 | STAR propensity | exact-X cells | $0.15\le e_x\le0.85$ | STAR 실행 금지 |
| G32 | STAR source law | frozen law | source-marginal error tolerance 충족 | STAR 실행 금지 |
| G33 | STAR role streams | each bundle | 역할 간 독립성 검증 | stage failure |
| G34 | STAR count | main | exactly 200 bundles, no extension | promotion 금지 |
| G35 | STAR labels | table·figure | 금지 metric·causal-truth 표현 없음 | publication 금지 |
| G36 | artifact integrity | each stage | manifest와 SHA256 전부 일치 | promotion 금지 |
| G37 | reproduction | tables·figures | raw에서 재계산한 값과 일치 | publication 금지 |
| G38 | claim map | manuscript handoff | 모든 문장에 허용 evidence mapping | 원고 수정 금지 |
| G39 | learner scope | CATE | DRF·RF, fixed spline $K_{\mathrm{spline}}=3$, neural 없음 | 실행 금지 |
| G40 | candidate grid | CATE | $\{0,.25,.5,.75,1\}^2$와 restriction subset 일치 | 실행 금지 |

Acceptance threshold가 이 표에 수치로 정의되지 않은 경우 구현자가 임의로 정하지
않는다. Numeric tolerance는 outcome을 보기 전 protocol fixture의 scale에 맞춰
고정하고 machine-readable config에 넣는다.

---

## 18. Current command와 planned command

현재 MM 실험에는 `mm_run.py`와 `mm_analyze.py`가 존재한다. 새 18-cell과 STAR
통합 명령은 아직 구현되지 않았다. 따라서 아래는 interface contract이지 현재
실행 가능한 명령이라는 주장이 아니다.

```text
# PLANNED: calculation-free preflight
python code/confirm/generalization_run.py --phase smoke --preflight

# PLANNED: synthetic phase
python code/confirm/generalization_run.py --phase smoke
python code/confirm/generalization_run.py --phase pilot
python code/confirm/generalization_run.py --phase main

# PLANNED: STAR phase
python code/confirm/generalization_star_run.py --phase smoke --star-csv <existing-file>
python code/confirm/generalization_star_run.py --phase pilot --star-csv <existing-file>
python code/confirm/generalization_star_run.py --phase main --star-csv <existing-file>

# PLANNED: immutable analysis only after promoted stages exist
python code/confirm/generalization_analyze.py --synthetic-stage <path> --star-stage <path>
```

Preflight는 scientific calculation과 결과 파일 생성을 하지 않고 다음만 반환한다.

- resolved protocol와 target stage
- code·protocol fingerprint
- 정확한 cell·seed schedule
- collision state
- required input과 output 목록
- STAR 사용 시 resolved path, hash와 schema status

실제 파일명과 CLI는 구현 PRD와 tests가 정할 수 있다. 변경하면 이 문서의 planned
interface와 대응표를 남겨야 한다.

---

## 19. 단계별 stop, resume, rollback

### 19.1 Smoke stop

다음 중 하나면 pilot으로 가지 않는다.

- cell 또는 role schedule 불일치
- ATE-CATE bundle fingerprint 불일치
- selection leakage
- ratio route가 다른 raw data를 사용
- mismatch perturbation이 $\tau(X)$를 변경
- hidden clipping, redraw 또는 failure omission
- STAR input·source-law preflight 실패

### 19.2 Pilot stop

다음 중 하나면 main으로 가지 않는다.

- 결과 schema가 raw에서 표·그림까지 연결되지 않음
- runtime 또는 memory 때문에 fixed main schedule을 완료할 수 없음
- failure가 명시적으로 분류되지 않음
- tail diagnostic이 extreme value를 잃거나 nonfinite를 숨김
- bootstrap이 method row를 독립 표본처럼 취급
- pilot 후 protocol을 바꿨지만 새 fingerprint와 새 main namespace를 만들지 않음

### 19.3 Main stop

Main 도중 hard validity failure가 발견되면 새 unit 시작을 중단한다. 이미 생성된
raw result와 failure ledger는 보존한다. 수정된 code로 이어 붙이지 않는다. 새
fingerprint와 새 stage에서 처음부터 실행한다.

### 19.4 Resume

같은 fingerprint에서만 resume한다. Ledger와 hash가 완전한 unit은 건너뛰고,
부분 unit은 같은 scheduled seed로 다시 시작한다. 다른 fingerprint의 result를
복사하거나 prefix로 재사용하지 않는다.

### 19.5 Rollback

Promotion 전 오류는 새 run이 소유한 scratch만 제거할 수 있다. Promoted stage,
기존 MM stage, 사용자 파일과 다른 scratch는 삭제하지 않는다. Promotion 뒤 오류는
기존 stage를 수정하지 않고 superseded marker를 별도 승인 아래 기록하며 새 stage를
만든다.

---

## 20. WHI/BioLINCC future feasibility

WHI는 synthetic과 STAR 결과 이후 고려할 수 있는 future natural-pair 후보다.
공식 자료 접근과 study documentation의 출발점은
[NHLBI BioLINCC WHI Clinical Trial and Observational Study](https://biolincc.nhlbi.nih.gov/studies/whi_ctos/)
이다.

이 PRD에서는 다음을 하지 않는다.

- 접근 가능하다고 가정
- treatment, outcome 또는 RCT-OBS pairing을 미리 확정
- STAR를 대신하는 confirmatory dataset으로 세기
- 접근 신청이나 다운로드를 현재 실행의 prerequisite로 만들기

WHI는 접근권, license·data-use 조건, 변수 dictionary, trial·observational cohort
mapping, estimand, support와 validation target을 별도 feasibility review가 승인한
뒤에만 새 protocol로 들어간다. 그러므로 현재 18-cell synthetic 또는 STAR
실행의 blocker가 아니다.

---

## 21. 최종 manuscript handoff

실행과 분석이 끝나면 원고 작성자에게 다음 packet만 넘긴다.

1. Promoted stage path와 SHA256 manifest
2. 정확한 독립 replication·failure count
3. Figure 1에서 Figure 3까지의 reproduction command와 source summary
4. Table 1과 Table 2의 machine-readable CSV
5. Claim-to-evidence matrix의 PASS·FAIL·N/A 상태
6. MM-main external consistency comparison
7. STAR conditional-resampling limitation 문구
8. theorem eligibility와 empirical-only label
9. protocol deviation과 unresolved issue 목록

원고 반영은 별도 작업이다. 다음 표현은 evidence가 있더라도 그대로 쓰지 않는다.

- “FusionPPI always improves risk.”
- “Estimated transport is cost-free.”
- “STAR verifies true heterogeneous treatment effects.”
- “Source-marginal equality establishes causal identification.”
- “ATE and CATE measure the same improvement.”

허용되는 문장은 실행한 SCM, $n_R$, regime, learner, ratio route와 uncertainty를
명시한다. 유리하지 않은 결과와 failure regime도 같은 가시성으로 보고한다.

---

## 22. 구현 순서

이 문서 승인 뒤의 별도 구현 작업은 다음 순서로 진행한다.

1. Machine-readable 18-cell protocol과 exact seed schedule 작성
2. Fixed-array ATE·CATE identity tests 작성
3. Unified role-bundle generator와 ATE-CATE fingerprint test 구현
4. Shared와 covariate-shift SCM2·SCM3 parity 구현
5. Outcome-mismatch DGP와 $h_k$, $\delta_k$, effect-invariance tests 구현
6. Ratio route pairing과 classifier diagnostics 구현
7. Synthetic preflight와 smoke 구현
8. STAR integrated bundle preflight와 smoke 구현
9. Raw·failure·resume·promotion infrastructure 검증
10. Analysis, 그림 3개, 표 2개의 fixture tests 구현
11. Smoke 실행과 독립 review
12. Pilot 실행과 runtime·schema review
13. 별도 confirmatory namespace에서 main 실행
14. Raw-to-summary reproduction과 claim audit
15. 별도 승인을 받은 경우에만 원고 반영

각 단계는 이전 단계의 hard acceptance를 통과한 뒤에만 시작한다. 실험의 효과
방향은 다음 단계로 진행할 engineering 조건이 아니다.

---

## 23. 완료 정의

이 연구 묶음은 다음이 모두 참일 때 완료다.

- 18-cell synthetic main의 정확히 3,600 cell-replications가 accounting됨
- STAR main의 정확히 200 conditional bundles가 accounting됨
- 운영 장부에 main source-data bundle 총 3,800개가 기록되되 synthetic과 STAR가
  추론에서 pooling되지 않음
- 모든 failure가 보존되고 redraw가 없음
- ATE와 CATE의 동일 bundle·독립 role 계약이 검증됨
- MM-main과 기존 ATE 결과가 primary inference에 섞이지 않음
- outcome mismatch와 transport route의 구현 identity가 통과함
- 필수 그림 3개와 표 2개가 immutable raw에서 재생성됨
- Engineering PASS와 scientific direction이 분리되어 보고됨
- STAR의 permitted·forbidden metric 경계가 지켜짐
- Claim-to-evidence audit가 완료됨
- WHI가 future feasibility로만 남음

이 정의는 FusionPPI가 반드시 우월해야 한다고 요구하지 않는다. 정확하게 설계된
실험이 개선, 무차이 또는 악화 중 무엇을 보이는지 손실 없이 기록하는 것이 완료다.

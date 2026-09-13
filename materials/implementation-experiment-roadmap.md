# DataFusionPPI 구현·실험 실행 로드맵

상태: 실행 로드맵 v1, 2026-09-11 기준.

이 문서는 연구 아이디어를 실제 코드, 단계별 실행, 검증 산출물, 논문 표와
그림으로 연결하는 실행 계획이다. 통계 공식과 고정 설계의 authority는
`materials/2026-09-09-implementation-handoff-JA.md` Version 3이다. 이 문서는
그 사양을 바꾸거나 다시 정의하지 않는다. 둘이 충돌하면 handoff v3가 우선하며,
충돌은 구현으로 우회하지 말고 protocol revision 대상으로 올린다.

## 1. 최종 목표와 완료의 의미

최종 과학 목표는 작은 무작위시험 RCT와 큰 관찰자료 OBS를 결합하는 두 통로가
언제 유효하고 언제 유용한지 분리해서 보이는 것이다.

- $\lambda$ 통로는 RCT AIPW pseudo-outcome 안의 outcome regression을
  RCT 적합치에서 OBS 적합치 쪽으로 이동시킨다.
- $\omega$ 통로는 OBS에서 얻은 treatment-effect prediction을 control variate로
  사용한다.
- 조건부 효과는 $\tau(x)=\mathbb E_R\{Y(1)-Y(0)\mid X=x\}$이고,
  ATE의 목표는 RCT 모집단 평균효과 $\theta=\mathbb E_R\{\tau(X)\}$이다.
- CATE candidate $\zeta(x)$는 $\tau(x)$를 추정하는 함수이며, 그 위험은
  $\mathcal R(\zeta)=\mathbb E_R[\{\zeta(X)-\tau(X)\}^2]$이다.

이하 $n_R$과 $N_O$는 study가 표시하는 RCT와 OBS sample budget이다. ATE에서는
각 source의 nuisance/tuning/evaluation으로 나누기 전 total이고, CATE에서는
study가 명시한 각 역할의 크기다. MSE는 평균제곱오차, RMSE는 그 제곱근,
MCSE는 반복실험 summary의 Monte Carlo standard error를 뜻한다.

완료는 단순히 11,000회가 끝났다는 뜻이 아니다. 다음이 모두 있어야 한다.

1. gate가 유효한 상태 토큰 하나를 기록한다.
2. `GO_FULL_TWO_CHANNEL`일 때만 전체 94-cell 실행이 열린다.
3. 각 예정 seed는 결과 또는 명시적 실패로 남고, 누락·교체·재추출이 없다.
4. 정리표와 그림이 replication-level artifact에서 재생된다.
5. 각 결과가 theorem check, diagnostic, empirical extension,
   real-data stress test 중 하나로 정확히 표시된다.
6. 부정적 결과와 혼합된 결과도 그대로 논문 근거로 남는다.

joint가 RCT-only보다 나쁘다는 사실 자체는 코드 오류가 아니다. 또한 이 연구는
추정된 계수의 보편적 no-harm을 주장하지 않는다.

## 2. Claim과 증거의 대응

| 질문 | 필요한 증거 | 허용되는 해석 |
| --- | --- | --- |
| 고정 $(\lambda,\omega)$ ATE 공식이 맞는가 | 지정 DGP의 평균·분산 identity와 oracle ordering | 가정이 맞는 cell의 theorem check |
| Algorithm 1이 유용한 계수를 찾는가 | RCT-only 대비 bias, variance, RMSE, coverage, coefficient recovery | finite-sample 성능과 calibration |
| covariate shift에서 ratio가 작동하는가 | drift, balance residual, ESS, normalization, coverage | identity check와 overlap diagnostic |
| CATE selector가 위험을 낮추는가 | select와 독립인 report 또는 독립 truth sample | 선택 후 성능 비교 |
| Theorem 6 bound가 성립하는가 | prespecified bounded SCM-B와 고정 $B_0,\bar r_0$ | 해당 bounded DGP에서만 theorem check |
| Theorem 7 quadratic가 보이는가 | spline, $\rho=0$, unclipped, $p<n_R^{\rm tune}$, condition 제한 | 조건을 만족한 row에서만 theorem check |
| STAR에서 pipeline이 coherent한가 | frozen source law, exact marginal, independent roles, report score | real-data stress test; causal truth가 아님 |

Gaussian mean shift와 STAR에서 Theorem 6 radius는 `N/A`다. Gaussian ratio는
common support를 가지지만 전역적으로 bounded하지 않는다. STAR에서
$P_O^X=P_R^X$와 $r_0=1$을 구성해도 causal identification이 자동으로 생기지 않는다.

## 3. 2026-09-11 현재 상태

### 3.1 완료된 G0

현재 canonical G0는 engineering·algebra smoke test까지만 완료했다.

| 항목 | 현재 증거 |
| --- | --- |
| protocol | `DataFusionPPI-v3-agile-gate-2026-09-10` |
| gate 상태 | `G0`, 결정 `INCONCLUSIVE` |
| performance 실행 | 8개 cell에서 replication ID 0 하나씩, 총 8회 |
| long-format 결과 | 372 metric rows |
| hard assertion | 32/32 pass |
| 전체 assertion | 39 pass, 4 `not_evaluable_at_g0` |
| test suite | 15/15 pass |
| 측정 runtime | 현재 구현과 장비에서 86.772687초 |
| code SHA-256 | `4d0cd7c87d3a23cf0b4449e5f69d94fcad5f43094f4805461a35698b96756e55` |
| handoff SHA-256 | `e1b2207beac31094f18856ed74391e1dc1f3bb665efbbf986f3f98231e950b39` |

372는 반복 횟수가 아니다. 하나의 replication에서 estimator·ratio·metric 조합별로
여러 행이 생겨서 얻은 long-format record 수다. 과학적 반복 단위는 8회다.

G0의 8개 performance cell은 다음과 같다.

- ATE-1 네 cell: reference, confounding 0/correct model,
  confounding 2/linear model, STAR real-$X$ reference.
- ATE-2 synthetic 두 cell: mild와 strong Gaussian shift.
- CATE-1 두 cell: SCM-B reference와 confounding 2, DRF spline $K=3$.

별도로 다음 두 smoke가 assertion만 만들었다.

- CATE-1 SCM-B reference의 neural 75-candidate path.
- STAR mathematics, $\alpha=0.8$, $n_R=200$의 CATE-2 source-law,
  DRF spline $K=3$, selection, reporting path.

이 둘은 performance replication이 아니고 재사용 prefix도 아니며 gate 600회 cap에
포함되지 않는다.

### 3.2 아직 검증되지 않은 것

현재 artifact의 다음 네 항목은 `not_evaluable_at_g0`다.

- population variance identity의 stochastic 3-MCSE 검사.
- exact-ratio cell의 oracle ordering 집계 검사.
- reporting score와 true risk difference의 stochastic 3-MCSE 검사.
- interval coverage 집계 검사.

현재 `code/v3_gate.py`는 G0만 실행한다. `--stage g1`은
`G1_NOT_IMPLEMENTED`로 fail closed하며, G2와 full-run entry point는 없다.
따라서 현재 산출물을 G1 또는 전체 실험의 증거라고 부르면 안 된다.

### 3.3 현재 실행 가능한 명령

다음은 현재 구현에서만 유효하다.

```bash
cd /Users/yonghanjung/paios/research/papers/DataFusionPPI/code
python3 -m unittest -v test_v3_gate.py
```

```bash
cd /Users/yonghanjung/paios/research/papers/DataFusionPPI
python3 code/v3_gate.py --stage g0 --artifact-dir materials --validate-artifacts
```

이미 보존된 G0를 이유 없이 다시 쓰지 않는다. 새 stage 구현은 별도 승인과
별도 artifact 보존 설계를 먼저 거친다.

## 4. 모든 stage가 지켜야 할 실행 계약

여기서는 실행자가 자주 어기는 경계만 요약한다. 정확한 식은 handoff v3 Section 3과
4를 그대로 구현한다.

| 영역 | 변경 불가 계약 |
| --- | --- |
| 역할 | ATE는 source별 nuisance/tuning/evaluation, CATE는 nuisance/tuning/selection/reporting을 독립 생성한다. 선택은 selection-only이고 frozen candidate 평가는 reporting-only다. |
| RNG | `design`, source-role별 namespace, `truth`, `oracle`, `gate_bootstrap`을 canonical tuple seed로 분리한다. |
| ATE score | $Z_0$은 RCT fit, $Z_1$과 양 source의 $\widehat g$는 OBS fit을 사용한다. 계수는 tuning에서 freeze하고 estimate·variance는 evaluation에서 계산한다. |
| oracle | 20k→40k→80k→최대 100k; $A,B$ tolerance는 $0.01s_0^2$, $C,D$는 $0.005s_0^2$. cap 실패는 MC-inconclusive이며 redraw하지 않는다. |
| ratio | normalization은 유한 diagnostic뿐이다. BAL 수렴은 residual $\le10^{-6}$이고 raw optimizer status와 분리한다. non-BAL balance는 `N/A`다. |
| Gaussian | $d^2=\log2,\log5$, $E_O(r_0^2)=e^{d^2}$, ESS $=e^{-d^2}$이며 sample maximum으로 bound를 만들지 않는다. |
| theorem scope | SCM-B clipped path만 고정 $B_0=9,\bar r_0=1$의 Theorem 6 check다. Theorem 7은 spline, $\rho=0$, unclipped, $p<n_R^{\rm tune}$, condition $\le10^8$만 해당한다. neural은 empirical extension이다. |
| CATE metric | synthetic은 독립 $10^5$ truth sample의 true risk, STAR는 independent reporting-score difference를 쓴다. |
| CATE-2 | exact $X$ tuple 전체와 $e_x\in[0.15,0.85]$ cell을 사용하고, OOF rank·weight를 사전 freeze하되 learner feature에서 제외한다. source-law equality와 realized share를 혼동하지 않는다. v2의 23-strata/rank-half 방식은 사용하지 않는다. |
| comparator | shrinkage/adaptive는 tuning-only이고 $C_{RO}$를 포함한다. naive pooling은 OBS pseudo-outcome을 직접 포함하며 heuristic label을 보존한다. |
| failure | 모든 seed는 결과 또는 explicit failure를 남긴다. 실패를 삭제·교체·redraw하지 않는다. |

## 5. Agile gate의 수치 장부와 결정 규칙

| panel | cell 수 | cell당 cap | cumulative maximum |
| --- | ---: | ---: | ---: |
| ATE-1 | 4 | 100 | 400 |
| ATE-2 synthetic | 2 | 50 | 100 |
| CATE-1 | 2 | 50 | 100 |
| 합계 | 8 |  | 600 |

G0는 cell당 1회이므로 누적 8회다. G1은 같은 8개 cell에서 ID 0--9,
cell당 10회, 누적 80회다. G0의 ID 0이 semantic prefix로 호환된다고 사전 판정된
경우에만 새로 72회를 실행한다. 호환되지 않으면 과거 G0는 보존하고 G1 정의에
따라 80회를 새로 실행한다. 현재 whole-file code hash 규칙만으로는 G1 구현 후의
재사용을 자동 승인할 수 없다.

Gate estimator set은 RCT-only, $\lambda$-only, $\omega$-only, joint와 정의된
경우의 oracle이다. Neural, NSW/PSID, full CATE-2와 full comparator set은 G1/G2
decision schedule에서 제외한다.

G2는 결정에 필요한 inconclusive cell만 incremental하게 늘린다. 모든 cell을 cap까지
자동으로 채우지 않는다. scientific relevant cell은 ATE-1 confounding 0/correct model과
CATE-1 SCM-B reference다.

각 relevant cell에서 ATE는 $L=$ MSE, CATE는 $L=$ true risk로 두고

$$
G_\omega=
\frac{\min\{\widehat L_{\mathrm{RCT}},\widehat L_{\lambda\text{-only}}\}
-\widehat L_{\mathrm{joint}}}
{\widehat L_{\mathrm{RCT}}}
$$

를 replication-level paired bootstrap으로 계산한다. 90% interval과
$\delta=0.05$를 사용한다. coverage guard는
$\operatorname{Coverage}_{\mathrm{RCT}}-
\operatorname{Coverage}_{\mathrm{joint}}$의 90% interval lower endpoint가
0.05를 넘으면 실패한다.

- hard validity/theorem/engineering failure: `REDESIGN`.
- hard failure가 없고, relevant cell 하나 이상에서 $G_\omega$ lower endpoint가
  0.05보다 크며 coverage guard 통과: `GO_FULL_TWO_CHANNEL`.
- GO가 아니고 relevant cell 중 하나라도 cap 미도달: `INCONCLUSIVE`.
- GO도 REDESIGN도 아니고 relevant cap 도달: `REFRAME`.

이 네 상태는 exhaustive하다. 10회짜리 G1은 방향, finite output, runtime을 보는
초기 단계일 뿐 coverage를 인증하거나 no-harm을 증명하지 않는다.

## 6. G1 결과를 보기 전에 확정할 사항

아래 항목은 outcome을 검사하기 전에 코드·fixture·memo에 freeze해야 한다.

1. paired bootstrap resample 횟수.
2. 90% interval의 정확한 convention: percentile, basic 등.
3. `gate_bootstrap` seed와 cell별 substep mapping.
4. coverage difference를 paired statistic으로 구성하는 정확한 row contract.
5. G0 prefix 재사용을 판단할 semantic implementation fingerprint.
6. G0, G1, G2 artifact를 동시에 보존할 versioning과 promotion 경로.

handoff는 bootstrap 원리, confidence level, namespace는 정하지만 resample 횟수와
interval convention은 정하지 않았다. 이를 결과를 본 뒤 선택하면 안 된다.

semantic fingerprint는 최소한 다음 구성요소를 포함해야 한다.

- protocol와 handoff hash.
- DGP, sample roles와 sizes, estimator, metric의 구현 fingerprint.
- seed canonicalization과 replication ID mapping.
- nuisance, ratio, oracle, CATE fit 함수의 semantic version.
- dependency pin과 데이터 source hash.

handoff의 기본은 code unchanged다. 따라서 어떤 code 변경도 보수적으로는 rep 0
재사용 불가이며 G1 80회를 새로 실행한다. 일부 함수 hash가 같다는 사실만으로 이
기본을 느슨하게 해서는 안 된다. G0 생성 code path의 동일성을 semantic fingerprint로
대체하려면 실제 G1 outcome 전에 prospective protocol revision과 명시적 승인을 받고,
fixture와 manifest로 판정을 검증해야 한다. 불명확하면 항상 rerun으로 fail closed한다.

현재 G0 네 artifact는 보존 대상이다. 미래 stage는 먼저 별도 temp 또는 stage-specific
경로에 완전하게 쓰고 검증한 뒤 publish한다. 새 실행 전에 기존 artifact를 삭제하는
방식은 금지한다. 최종 경로 이름은 G1 구현 승인 전에 정할 미해결 결정이다.

## 7. 구현 backlog: 작은 end-to-end milestone

| milestone | 최소 산출물 | observable acceptance |
| --- | --- | --- |
| M0 G0 보존·compatibility | immutable G0 snapshot, semantic fingerprint, stage-specific atomic promotion | 실패한 새 실행이 G0를 바꾸지 않고, DGP·seed·metric 변경이 incompatibility로 검출되며 reuse 근거가 machine-readable함 |
| M1 decision 선구현 | frozen bootstrap config, paired $G_\omega$, coverage difference, exhaustive state machine, fixture tests | 실제 G1 outcome 전에 freeze; 수작업 fixture 일치; 0.05 경계의 strict inequality와 hard-failure 우선순위 확인 |
| M2 G1 engine | ID 0--9 executor, exact schedule validator, finite/runtime/role checks | invalid G0에서 fail closed; 80회 모두 result 또는 explicit failure; neural·NSW/PSID·full CATE-2·full comparator 제외 |
| M3 실제 G1 | cumulative 80회와 최초 gate memo | 승인된 reuse면 72회 추가, 아니면 과거 G0를 보존하고 80회 새 실행; 10회로 coverage/no-harm 주장 금지 |
| M4 targeted G2 | inconclusive-cell manifest와 incremental results | ATE-1 cap 100, ATE-2/CATE-1 cap 50; prefix ID 유지; REDESIGN이면 성능 실행 중단 |
| M5 full lock·runner | GO-only entry point, resumable manifest, seed ledger, per-study atomic artifact | 다른 세 token은 fail closed; resume 중복 없음; 실패 seed 보존 |
| M6 ATE-1 | 44 cells/4,400회와 네 baseline, reference-only crossfit panel | identity·oracle·Algorithm 1·interval을 분리하고 comparator tuning freeze/$C_{RO}$ 검사 |
| M7 sentinels | ATE-1 직후 네 reference cell에 각 400회 추가 | 총 1,600회; IDs 100--499; $[0.921,0.979]$ 밖이면 diagnose/no-go, redraw 금지; unresolved no-go 뒤 후속 study 금지 |
| M8 ATE-2 synthetic | 12 cells/1,200회와 다섯 ratio | analytic ratio moments·ESS와 empirical diagnostic 병기; BAL tolerance; clipping drift label; global-bound 문구 없음 |
| M9 CATE-1 | 24 cells/2,400회, spline 다음 neural | 6 learner/sieve 조합은 cell 수에 미가산; select/report 독립; Theorem 6/7/neural label 분리 |
| M10 NSW/PSID | 6 cells/600회, CPS와 PSID 모두 | URL/hash/row count 검사; point/length/drift/ESS; truth RMSE나 oracle ratio를 꾸미지 않음 |
| M11 CATE-2 | 8 cells/800회 | outcome-$\alpha$별 source law 사전 freeze; table fingerprints; STAR report score 사용; causal claim 금지 |
| M12 paper package | frozen CSV 기반 summary, figures, tables, study memos | MCSE와 label 유지; unstable denominator는 paired difference; negative/mixed finding 보존 |

세부 순서는 중요하다. M1의 statistic·decision 코드와 fixture를 먼저 review하고,
그 뒤 M2를 완성한 다음에만 M3에서 실제 G1 outcome을 연다. 결정성 검사는 작은
fixed-seed fixture 또는 한 cell로 수행하며 전체 80회 batch 두 번을 자동 요구하지 않는다.

CATE-1의 6개 조합은 DRF/RF $\times$ spline $K=3$, spline $K=5$, neural이다.
Spline은 25개 $(\lambda,\omega)$ 후보, neural joint는 75개
$(\lambda,\omega,\rho)$ 후보를 사용하지만 모두 기존 24개 design cell 안의 fit이다.

NSW의 참조값 1794.3은 full experiment의 descriptive point difference이지 알려진
truth가 아니다. CATE-2는 mathematics와 reading을 분리하고 true CATE risk 대신
independent reporting-score difference를 쓴다.

M7의 $[0.921,0.979]$ 판정은 nominal 0.95, 500회라는 계산과 함께
exact-ratio·correctly-specified sentinel의 coverage acceptance에 적용한다. STAR처럼
empirical-extension인 row는 같은 수치를 보고하더라도 theorem assertion으로 바꾸지
말고 handoff의 분류에 따라 calibration 또는 stress diagnostic으로 해석한다.

## 8. 전체 실행 수 장부

| study | 고정 design과 sample | unique cells | screening replications |
| --- | --- | ---: | ---: |
| ATE-1 | SCM 1/2/3/STAR $\times$ 11; reference는 $n_R=100,N_O=5000$, confounding 1, correct model. $N_O=5000$에서 $n_R\in\{50,100,200,400\}$, $n_R=100$에서 $N_O\in\{1000,20000\}$ 추가, quality는 $c\in\{0,1,2\}\times\{\text{correct,linear}\}$에서 중복 reference $(1,\text{correct})$를 뺀 5개 추가 | 44 | 4,400 |
| ATE-2 synthetic | 3 SCM $\times$ $n_R=100,400$ $\times$ $d^2=\log2,\log5$; $N_O=5000$ | 12 | 1,200 |
| ATE-2 NSW/PSID | NSW treated 185명을 반으로 나누고 RCT 쪽에 NSW controls 260명을 결합; $n_R=100,200,$ all $\times$ CPS/PSID | 6 | 600 |
| CATE-1 | SCM-B/SCM 2/SCM 3/STAR $\times$ 6; reference는 $n_R=200,N_O=5000$, confounding 1, base $\tau$. 축은 $n_R\in\{100,200,400\}$, reference 주위 quality 추가 $c\in\{0,2\}$, strong variant $\tau\mapsto2\tau$ 추가 1개 | 24 | 2,400 |
| CATE-2 | $n_R=200,400$ $\times$ $\alpha=0.4,0.8$ $\times$ mathematics/reading | 8 | 800 |
| performance 합계 |  | 94 | 9,400 |
| ATE-1 sentinel 추가 | 기존 4 reference cells에 각 400회 | 기존 cell | 1,600 additional |
| unique full-plan 합계 |  | 94 | 11,000 |

ATE의 기본 source별 role 비율은 RCT nuisance/tuning/evaluation
$0.4/0.3/0.3$, OBS $0.6/0.2/0.2$이며 study별 명시가 우선한다. CATE-1의 표시된
$n_R,N_O$는 nuisance, tuning, selection, reporting 각각의 독립 sample 크기다.
CATE-2는 nuisance와 tuning에 표시된 $n_R$와 $N_O=5000$을 각각 쓰고,
selection과 reporting은 source별 1,000개씩 독립 추출한다.

ATE-1은 RCT-only, $\lambda$-only, $\omega$-only, joint,
oracle-coefficient joint의 다섯 유형에 semi-supervised, shrinkage, adaptive,
naive pooling을 비교한다. ATE-2 synthetic은 다섯 ratio를 비교하고 real panel은
classifier, BAL-X, BAL-X+$g$만 쓴다. CATE-1/2는 RCT-only, $\lambda$-only, $\omega$-only,
joint, grid oracle 및 DRF oracle calculus와 feasible plug-in, experimental grounding,
domain-indicator pooling을 같은 역할과 sieve에서 비교한다.

11,000은 unique design replications다. estimator 수, ratio 수, 25/75 candidate,
6개 learner/sieve 조합, long-format metric row 수, model fit 수가 아니다.

gate의 최대 600회는 full-study prefix로 실제 재사용된 부분에 한해서만 중복 계산하지
않는다. protocol 밖 stress corner는 별도 non-reusable 실행으로 표시하고 더한다.

handoff의 9--10 single-core hour는 잠정 예산이다. 현재 G0 runtime을 11,000회에
단순 선형 외삽하지 않는다. G1에서 역할·learner별 시간을 측정하고, parallel scaling을
직접 측정하기 전에는 8-core 소요시간을 약속하지 않는다.

## 9. 계획된 산출물

### 9.1 Gate

현재 G0 artifact:

- `materials/gate_v3_replications.csv`
- `materials/gate_v3_summary.csv`
- `materials/gate_v3_assertions.json`
- `materials/gate_v3_results.md`

G1/G2 stage-specific 보존 경로는 M0에서 정한다. 현재 G0를 덮어쓰는 이름을
그대로 사용하지 않는다.

### 9.2 Study artifacts

- ATE-1: `ate1_when_fusion_helps_replications.csv`, summary CSV,
  main figure, coefficient-recovery figure.
- ATE-2 synthetic: `ate2_shift_replications.csv`, summary CSV, figure.
- ATE-2 real: `ate2_nsw_replications.csv`, summary CSV, figure.
- CATE-1: `cate1_sieve_validation_replications.csv`, summary CSV,
  validation, regret, coefficient figures.
- CATE-2: `cate2_star_real_replications.csv`, summary CSV, figure.
- 각 study에는 design, 실제 실행 수, measured compute, 실패, 분류 label,
  3--5문장 결과를 담은 Markdown memo를 붙인다.

### 9.3 논문 본문과 부록

| 위치 | 내용 |
| --- | --- |
| Figure 1 | ATE-1의 RCT-only 대비 RMSE, correct/misspecified 두 panel |
| Figure 2 | CATE-1 risk, DRF/RF와 spline/neural |
| Table 1 | ATE coverage·length와 bounded-DGP CATE regret |
| Table 2 | NSW와 STAR real-outcome 두 panel |
| Appendix | ratio 비교, coefficient recovery, exact variance, leading risk, honest-vs-crossfit |

그림은 one-column width, 약 55 mm height, 최종 크기에서 읽히는 font로 만든다.

## 10. Stage별 stop 조건과 복구 원칙

즉시 중단 조건:

- handoff, data, code, dependency fingerprint가 manifest와 다름.
- role aliasing, selection/report leakage, global RNG 사용.
- exact identity나 source-law hard assertion 실패.
- 예정 schedule cardinality와 실제 row가 다름.
- 결과가 유한하지 않거나 seed가 누락됨.
- unauthorized protocol change가 발견됨.

분류해서 처리할 조건:

- BAL convex-hull infeasibility: explicit infeasible, 삭제·redraw 금지.
- oracle MCSE cap 실패: Monte Carlo-inconclusive, 같은 replication ID 유지.
- CATE-2 directional bias diagnostic 실패: soft diagnostic; hard engineering failure 아님.
- joint underperformance: scientific result; 구현 실패 아님.
- sentinel coverage 이탈: diagnosis/no-go signal; redraw 금지.

복구는 마지막으로 검증된 stage artifact에서 시작한다. 새 실행은 별도 경로에서
완성·검증한 뒤 publish하며, 이전 성공 artifact를 먼저 삭제하지 않는다.

## 11. Planned CLI contract

아래는 구현 목표이며 현재 실행 가능한 명령이 아니다.

```bash
# PLANNED ONLY — G1 implementation과 frozen config가 먼저 필요
python3 code/v3_gate.py --stage g1 --artifact-dir <stage-specific-dir> \
  --gate-config <frozen-gate-config>
```

```bash
# PLANNED ONLY — decision-directed extension
python3 code/v3_gate.py --stage g2 --cells <inconclusive-cell-manifest> \
  --artifact-dir <stage-specific-dir>
```

```bash
# PLANNED ONLY — GO token과 validated manifest 없이는 fail closed
python3 code/run_full_v3.py --decision-memo <validated-gate-memo> \
  --schedule <full-schedule-manifest> --resume
```

실제 option 이름과 새 파일 생성은 각 구현 turn의 승인 후 확정한다. roadmap에
명령 예시가 있다는 사실은 코드나 새 artifact 생성 승인이 아니다.

## 12. 최종 QA checklist

- [ ] G0 artifact와 manifest가 immutable하게 보존됨.
- [ ] bootstrap 횟수·interval convention·seed가 outcome 전에 freeze됨.
- [ ] semantic fingerprint reuse fixture가 통과함.
- [ ] G1 누적 80회와 exact schedule이 검증됨.
- [ ] G2가 필요한 cell만 늘리고 600 cap을 지킴.
- [ ] decision token이 정확히 하나임.
- [ ] GO가 아니면 full run이 fail closed함.
- [ ] 94 cells, 9,400 screening, 1,600 extension, 11,000 total이 일치함.
- [ ] 네 sentinel 모두 500회이며 interval 기준은 $[0.921,0.979]$임.
- [ ] ATE true $\theta$, synthetic CATE true risk, STAR report score가 구분됨.
- [ ] CPS와 PSID가 모두 포함됨.
- [ ] CATE의 6 learner/sieve 조합을 unique cell로 중복 계산하지 않음.
- [ ] reference-only honest-vs-crossfit panel이 포함됨.
- [ ] normalization diagnostic과 BAL convergence가 분리됨.
- [ ] Gaussian과 STAR에 sample-max global bound 문구가 없음.
- [ ] neural result가 Theorem 7 표에 들어가지 않음.
- [ ] CATE-2 source equality를 causal identification으로 부르지 않음.
- [ ] 모든 seed가 result 또는 explicit failure를 가짐.
- [ ] runtime과 parallel scaling은 측정값과 미측정 추정을 구분함.
- [ ] figure와 table가 replication artifact에서 재생됨.

## 13. 아직 결정해야 할 항목

다음은 현재 미해결이며 결과를 보기 전에 별도 승인으로 정한다.

1. gate paired bootstrap의 resample 수와 90% interval convention.
2. semantic compatibility fingerprint의 exact field와 승인 절차.
3. G0/G1/G2 immutable artifact의 최종 파일명과 promotion 방식.
4. full executor의 checkpoint granularity와 parallel backend.
5. G1 측정 뒤 갱신할 compute budget과 worker 수.
6. baseline 비교 목록은 handoff대로 고정하고, 비차단 pending 항목인
   predictions-as-surrogates related-work citation을 재확인할 시점.

이 항목을 임의 threshold나 사후 선택으로 메우지 않는다. 결정이 필요한 시점에는
선택지, 결과에 미치는 영향, handoff와의 호환성을 함께 제시한다.

## 14. 권장 실행 순서 요약

1. 현재 G0를 보존하고 semantic compatibility 계약을 만든다.
2. bootstrap과 decision engine을 fixture로 먼저 구현·freeze한다.
3. G1 executor와 validator를 구현한 뒤 누적 80회를 실행한다.
4. 필요할 때만 G2를 늘려 terminal decision을 기록한다.
5. GO일 때만 full executor를 연다.
6. ATE-1 screening 뒤 네 sentinel을 500회로 확장하고 no-go를 먼저 해소한다.
7. 그 뒤 ATE-2 synthetic, CATE-1 spline, CATE-1 neural,
   ATE-2 NSW/CPS·PSID, CATE-2 순으로 실행한다.
8. frozen CSV에서 summary, figures, tables, study memo를 재생한다.
9. claim과 evidence label을 점검한 뒤에만 manuscript 반영을 별도 승인받는다.

이 순서는 계산량보다 먼저 식별 가능성, honest evaluation, 재현성, 실패 보존을
지킨다. 이번 작업은 roadmap 한 파일만 만든다. 향후 승인된 구현 범위 안의 routine
step은 이 순서를 따르되, 통계 protocol deviation과 commit·push는 별도 지시가 필요하다.

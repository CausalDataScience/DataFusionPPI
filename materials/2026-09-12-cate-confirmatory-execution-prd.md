# DataFusionPPI CATE 원고 충실 최소 실험 PRD

개정일: 2026-09-13. 상태: **계획 업데이트 완료, 구현·실험 미수행**.

## 1. 목표와 적용 계약

최종 목표는 원고의 DRF와 RF CATE 추정기를 구현하고, 공변량 분포가 같은 경우와 다른 경우에 같은 RCT 목표의 위험·편향·분산을 비교하는 것이다. 용한과 후속 구현자는 이 문서로 구현과 실험의 완료를 판정한다. 산출물은 검증된 공통 엔진, 두 실험 셀의 원시 결과, 그림 2개, 결과표 1개, 재현 명령과 설정이다.

통계적 authority는 현재 `manuscript/main/2.tex`, `3.tex`, `4.tex`다. 이 PRD가 기존 handoff 및 아래 접힌 구 실행계획의 실행 범위를 대체한다. 구 M0–M5 승격, 35개 검사, STAR 재구성, 네 사분면 보정은 현재 실행의 선행조건이 아니다. 원고 `main/5.tex`의 넓은 CATE roadmap은 후속 확장으로 남긴다. 현재 범위 완료를 원고 전체 실험 또는 학회 제출 완료라고 부르지 않는다.

이번 요청에서 변경하는 것은 계획 문서다. 코드, 실험 실행, 원고, mirror, commit, push는 이번에 수행하지 않는다. 후속 구현은 기존 canonical project 파일을 최대한 재사용한다.

## 2. 계속 적용할 필요성 검사

구현 전, smoke 후, pilot 후, 본 실험 후 다음 질문과 답을 짧게 남긴다.

1. 이 작업이 원고의 어느 식 또는 아래 실험 질문에 필요한가?
2. 기존 함수·자료·후보 예측을 재사용하여 답할 수 있는가?
3. 빼면 구현 정확성이나 결과 해석이 불가능해지는가?

세 번째 답이 아니오이면 이번 범위에서 제외한다. 좋은 결과를 얻기 위한 DGP, seed, 교란, 기저 탐색은 하지 않는다. 이득이 없다는 결과도 유효한 완료다.

| 필요한 질문 | 최소 증거 | 제외하는 추가 작업 |
| --- | --- | --- |
| 공통 법칙에서 두 채널이 무엇을 바꾸는가? | shared 1셀, 네 후보 부분집합 | 네 사분면·손잡이 보정 |
| 정확한 전송비율을 쓰면 어떤 유한표본 성능을 얻는가? | shifted-exact 경로 | 여러 이동 강도 |
| 추정 가능한 비율을 쓰면 성능이 얼마나 달라지는가? | 같은 shifted 자료의 classifier 경로 | ratio 추정기 대회 |
| 이동을 무시하면 어떤 문제가 생기는가? | 같은 자료의 wrong-unit 경로 | 새 실패 DGP |
| 위험 변화가 분산 때문인가 편향 때문인가? | 같은 예측의 반복 분해 | 별도 분산 실험 |

## 3. 원고–구현 대응

$X$는 처치 전 공변량, $A$는 처치 여부, $Y$는 결과다. RCT는 무작위시험, OBS는 관측자료다. 목표 CATE는 $\tau(x)\triangleq\mathbb E_R\{Y(1)-Y(0)\mid X=x\}$, 위험은 $\mathcal R(\zeta)\triangleq\mathbb E_R(\zeta-\tau)^2$다.

| 기능 | 원고의 정확한 locator | 기존 구현 |
| --- | --- | --- |
| 자료 분리·목표 | `2.tex`: `ass:honest-separation`, `eq:cate` | 역할 분리와 SCM |
| 혼합 회귀·의사결과 | `2.tex`: `def:fusion-score`; `3.tex`: `eq:obs-effect-prediction` | `build_scores`, `ghat_on` |
| 두 적합 손실 | `4.tex`: `eq:dr-fusion-loss`, `eq:rf-fusion-loss` | `learner_parts`, `sieve_solver` |
| 후보 선택 | `4.tex`: `alg:cate-fusion-selection`, `eq:common-reference-validation` | `validation_score` |
| 선형 기저 해 | `4.tex`: `eq:sieve-learner`, `eq:sieve-normal-equation`, `eq:sieve-lambda-path` | 공통 선형 solver |
| 정확·추정 비율의 차이 | `4.tex`: `prop:cate-target`, `eq:common-reference-identity` | ratio 경로와 진단 |
| 분류기 비율 | `3.tex`: `eq:ratio-classifier-estimator` | prior-corrected classifier |

### 3.1 적합식

$\widehat\mu_R,\widehat\mu_O$는 nuisance 표본에서 학습한 결과 회귀다. 기존 `fit_outcome(spec="flexible", alpha=5.0)`를 두 원천에 동일하게 적용한다. $V=(X,A,Y)$와 회귀쌍 $q=(q_0,q_1)$에 대해 원고 `def:fusion-score`의 $\varphi(V;q)\triangleq q_1(X)-q_0(X)+A\{Y-q_1(X)\}/e(X)-(1-A)\{Y-q_0(X)\}/\{1-e(X)\}$를 쓴다. $Z_0=\varphi(V;\widehat\mu_R)$, $Z_1=\varphi(V;\widehat\mu_O)$, $Z_\lambda=(1-\lambda)Z_0+\lambda Z_1$로 둔다. $\widehat g=\widehat\mu_{O,1}-\widehat\mu_{O,0}$다. RF에서는 알려진 RCT 성향점수 $e$를 써서 $m_S=e\widehat\mu_{S,1}+(1-e)\widehat\mu_{S,0}$, $m_\lambda=(1-\lambda)m_R+\lambda m_O$, $\chi=(A-e)^2/\{e(1-e)\}$로 둔다.

다음 경험평균은 tuning 표본에서 계산한다.

$$
\widehat L_{\mathrm{DRF}}(\zeta)=\mathbb P_R(Z_\lambda-\zeta)^2+\omega\{\mathbb P_O[r(\widehat g-\zeta)^2]-\mathbb P_R(\widehat g-\zeta)^2\}.
$$

$$
\widehat L_{\mathrm{RF}}(\zeta)=\mathbb P_R\frac{[Y-m_\lambda-(A-e)\zeta]^2}{e(1-e)}+\omega\{\mathbb P_O[r(\widehat g-\zeta)^2]-\mathbb P_R[\chi(\widehat g-\zeta)^2]\}.
$$

두 식에서 $\zeta(x)=b(x)^\top\beta$로 두고 $\rho\|\beta\|^2$를 더해 최소화한다. RF 적합의 RCT 보정항에는 $\chi$가 반드시 있다. 기존 spline3 한 가지와 고정 $\rho=10^{-2}$를 모든 learner·경로에 적용한다. Knots와 표준화는 RCT nuisance만 사용해 tuning 전에 고정한다. 반복마다 nuisance가 달라지므로 기저도 다시 학습한다. 결과는 회귀·기저·비율·후보 선택을 포함한 전체 절차의 성능이다.

정상 행렬에는 선언한 ridge 외 jitter를 추가하지 않는다. 적합·선택·평가에서 $\widehat\zeta$, $\widehat g$, $Z_0$를 임의 clipping하지 않는다. 해가 비유한이거나 정규방정식 잔차를 만족하지 않으면 실패로 기록한다. Gaussian 결과와 비율은 전역적으로 유계가 아니므로 정리 6 Hoeffding 반경은 `N/A`다. 정리 7의 대수 검산은 가능하지만 $\rho>0$ 결과로 $\rho=0$ 또는 점근 remainder를 검증했다고 주장하지 않는다.

### 3.2 후보 선택

공통 격자는 $\lambda,\omega\in\{0,0.5,1\}$의 9개다. 각 learner·ratio에서 한 번 적합한 후보를 다음 네 부분집합으로 읽는다.

| 규칙 | 후보 부분집합 |
| --- | --- |
| RCT-only | $(0,0)$ |
| $\lambda$-only | $\omega=0$ |
| $\omega$-only | $\lambda=0$ |
| joint | 전체 9개 |

모든 규칙은 독립 selection 표본의 원고 점수를 최소화한다. $G=(\widehat g-\widehat\zeta)^2$는 각 후보에서 계산한다.

$$
\widehat C_j(\lambda,\omega)=\mathbb P_R^{\mathrm{sel}}[(Z_0-\widehat\zeta)^2-\omega G]+\mathbb P_O^{\mathrm{sel}}[\omega rG].
$$

선택의 $\omega$는 적합 후보의 $\omega$와 같다. 별도 validation-$\omega$는 없다. RF도 이 점수를 쓰며 selection RCT 보정항에는 $\chi$를 넣지 않는다. 동점은 $(\lambda,\omega)$ 사전식 순서로 정한다. 모든 부분집합이 같은 ridge·자료를 쓴다. 참 $\tau$와 truth risk를 학습·선택에 전달하지 않는다. 평가 시 같은 9개 중 true risk 최소 후보를 `grid_oracle`로만 표시한다. 별도 oracle 최적화는 하지 않는다.

## 4. 고정 실험 설정: 두 셀만

`fusion_data.scm_family(1, confounding=1, heterogeneity=1)`의 구조 계수를 고정한다. SCM1 하나만 사용한다. 기존 공분산 $\Sigma$, 원래 평균 $\mu$, 이동 벡터 $v$를 읽고 $d\triangleq v/\sqrt{v^\top\Sigma^{-1}v}$로 정한다. 이동의 Mahalanobis 크기는 1이다. $v=0$이면 탐색 없이 $v=(1,0,0,0)^\top$로 대체한다.

| 셀 | RCT 공변량 | OBS 공변량 | ratio 경로 |
| --- | --- | --- | --- |
| shared | $N(\mu,\Sigma)$ | $N(\mu,\Sigma)$ | 정확한 $r_0=1$ |
| shifted | $N(\mu,\Sigma)$ | $N(\mu+d,\Sigma)$ | exact / classifier / wrong-unit |

Shifted exact 비율은

$$
r_0(x)=\exp\{-d^\top\Sigma^{-1}(x-\mu)+\tfrac12d^\top\Sigma^{-1}d\}.
$$

두 셀 모두 알려진 RCT 성향점수를 **$e=0.35$**로 정한다. $e=0.5$에서는 DRF와 RF가 관측값별로 같아지므로, 실제 두 추정기를 구분하는 최소 설정이다. $\tau$, 잠재 교란 구조, outcome 식과 잡음은 기존 SCM과 같다. 기존 생성식을 재사용하면서 OBS 평균만 바꾸고, 바뀐 $X$에서 OBS propensity·처치·결과를 다시 생성한다. RCT 역시 $e=0.35$로 처치를 생성한 뒤 기존 outcome 식을 적용한다. 기존 $Y$를 둔 채 $X$만 이동시키지 않는다. RCT target을 이동시키는 기존 `exp_cate3.py`를 그대로 비교 driver로 쓰지 않는다. 추가 propensity trimming은 하지 않는다.

| source | nuisance | tuning | selection | 총수 |
| --- | ---: | ---: | ---: | ---: |
| RCT | 100 | 100 | 100 | 300 |
| OBS | 3,000 | 1,000 | 1,000 | 5,000 |

각 source의 세 역할은 독립이며 겹치지 않는다. Synthetic reporting은 실행 전 별도로 생성·고정한 RCT 목표의 **50,000개 truth-grid 공변량과 참 $\tau$**로 수행한다. 모든 평균 위험·편향·분산은 정확히 이 grid에서 계산한다. 추가 labeled reporting 표본은 필요 없다.

Classifier는 기존 standardized logistic regression 한 가지다. nuisance 학습 prior $\widehat\pi=100/(100+3000)$로 $\widehat r=(1-\widehat\pi)\widehat s/\{\widehat\pi(1-\widehat s)\}$를 쓴다. `decision_function`의 log odds를 이용해 probability clipping 없이 계산한다. 기존 `ratio_classifier`의 숨은 probability clipping은 재사용 전에 수정한다. Ratio clipping, selection 표본 정규화, oracle로 추정 비율 보정을 하지 않는다. overflow나 미수렴은 명시적 실패다.

반복 ID 하나에서 두 셀의 RCT 배열을 공유한다. 같은 OBS 잠재변수와 기본 오차에 평균 이동을 적용하여 셀 간 비교도 pairing한다. 각 셀 내 모든 ratio·learner는 동일 자료·회귀·기저·selection을 공유한다. Route 이름을 자료 RNG seed에 넣지 않는다. $\omega=0$ 예측은 경로 간 같아야 한다. 서로 다른 반복 ID만 독립이다.

## 5. 재사용과 필요한 수정

| 기존 코드 | 작업 |
| --- | --- |
| `fusion_data.py`, `ssem_ate_pilot.py` | SCM과 결과식을 재사용; OBS 평균 및 공통 $e=0.35$의 작은 생성 어댑터 |
| `fusion_core.py` | `build_scores`, `learner_parts`, `sieve_solver`, `validation_score`를 단일 수식 구현으로 사용; 숨은 jitter·probability clipping 수정 |
| `fusion_cate.py` | 두 learner·부분집합 선택 재사용; clipping, 고정 split, 외부 비교군·불필요 oracle 호출 비활성화 또는 수정 |
| `confirm/cate_engine.py`, `cate_run.py`, `cate_protocol.py` | 역할·실행 관리 재사용; 복제 대수·DRF 전용 solver·RCT-only 선택·true-risk ablation 선택을 공통 함수 호출로 교체 |
| `fusion_io.py`, 기존 집계·plot | 아래 정의와 맞는 저장·분해·재표집만 재사용 |
| 기존 테스트·M1 fixture | 수식 검산만 재사용; 전체 승격 체계 재실행은 선행조건에서 제외 |

새 프레임워크나 새 대수 엔진은 만들지 않는다. 기존 산출물은 식·역할·목표법칙·code provenance가 이 PRD와 같음을 확인한 때만 재사용한다. 좋은 과거 수치를 새 main에 합치지 않는다. ATE는 재실행하지 않는다.

## 6. 구현과 실행 순서

다음은 **향후 기존 runner에 구현할 CLI 계약**이다. 현재 실행 가능하다는 뜻이 아니다.

```bash
cd /Users/yonghanjung/paios/research/papers/DataFusionPPI/code
python3 -m confirm.cate_run --plan manuscript-minimal-v1 --phase preflight
python3 -m unittest -v test_cate_confirm
python3 -m confirm.cate_run --plan manuscript-minimal-v1 --phase smoke
python3 -m confirm.cate_run --plan manuscript-minimal-v1 --phase pilot
python3 -m confirm.cate_run --plan manuscript-minimal-v1 --phase main
```

| 단계 | 셀당 반복 | 두 셀 record 수 | 목적 |
| --- | ---: | ---: | --- |
| 구현·고정 배열 검사 | 0 | 0 | §8의 원고 일치 검산 |
| smoke | 1 | 2 | 두 learner와 모든 경로가 끝까지 연결 |
| pilot | 20 | 40 | 오류·runtime·꼬리 확인, 본 실행 feasibility |
| main | 200 | 400 | 고정 예산의 최종 성능 비교 |

Main은 **200개 독립 seed 묶음, 400개 셀·반복 record**다. 두 셀 RCT를 공유하므로 400회 완전 독립 실험이라고 부르지 않는다. Learner×route 묶음은 $(1+3)\times2\times200=1,600$개, 네 규칙 결과는 6,400행, 9개 격자의 논리적 후보 적합 수는 14,400개다. $\lambda$별 동일 행렬과 중복 $\omega=0$ 후보 재사용으로 실제 선형 계산은 줄인다. 행 수를 반복 수로 세지 않는다.

Seed는 protocol/phase/replication/source/role로 정하고 phase별로 분리한다. Smoke와 pilot은 main에 합치지 않는다. Pilot은 구현 오류와 계산 가능성을 판단하며 좋은 방향이 나와야 main으로 가는 gate를 두지 않는다. DGP·모형·격자를 성능에 맞춰 변경하지 않는다.

Main은 셀당 200회에서 종료한다. 평균 차이의 Monte Carlo 표준오차(MCSE)는 paired 반복차이의 표준편차를 $\sqrt{200}$으로 나눈다. 정밀도가 부족하면 `INCONCLUSIVE`와 MCSE를 보고한다. 자동 500회 확대, 유의성이 나올 때까지 반복, 결과 방향에 따른 중단은 없다. 후속 반복은 해결할 불확실성과 기존 결과만으로 답할 수 없는 이유를 제시한 별도 요청으로 다룬다.

## 7. 위험·편향·분산과 최소 산출물

반복 수를 $K$, 반복별 예측을 $f_k(x)$, 반복 평균을 $\bar f(x)$라 하자. $\mathbb P_T$는 고정 truth grid의 평균이다.

$$
\overline R=K^{-1}\sum_k\mathbb P_T(f_k-\tau)^2,\qquad
\widehat B^2=\mathbb P_T(\bar f-\tau)^2,\qquad
\widehat V=(K-1)^{-1}\sum_k\mathbb P_T(f_k-\bar f)^2.
$$

$$
\overline R=\widehat B^2+\frac{K-1}{K}\widehat V.
$$

분산은 같은 $x$에서 독립 반복별 예측의 변동이다. 위험 수치의 분산이나 서로 다른 $x$ 사이 예측값 분산을 대신 쓰지 않는다. $\widehat B^2$는 유한 반복의 plug-in 제곱편향이다.

Primary는 joint 위험에서 같은 learner RCT-only 위험을 뺀 paired 평균이다. Joint 대 $\lambda$-only 차이도 보고하여 $\omega$의 추가 효과를 읽는다. 네 규칙 모두 평균 위험·분산·제곱편향·baseline 대비 차이와 비율을 낸다. 위험비는 **평균 위험의 비율**, 분산비는 위 $\widehat V$의 비율이다. 분모 0은 `N/A`다. 중앙값으로 평균을 대체하지 않는다. 90/99분위수와 최댓값은 꼬리 진단이다.

95% 구간은 seed 묶음을 통째로 재표집하는 paired bootstrap 2,000회, seed 91307로 계산한다. 두 셀·모든 경로·learner·truth 예측을 같은 재표집에 묶는다. 점추정과 구간이 같은 집계 함수를 사용한다. 구간은 유한 truth grid에 조건부인 근사 Monte Carlo 구간이며 여러 비교의 동시 95% 보장을 주장하지 않는다. 실패를 0이나 새 seed로 대체하지 않는다. 비교는 양쪽이 모두 성공한 반복 ID의 교집합만 사용하고 paired 수와 각 방법의 전체 실패율을 함께 보고한다. 분해 항등식과 해당 비교의 위험·편향·분산은 동일한 성공 ID 집합에서 계산한다. 성공 반복 조건부 요약을 무조건부 성능으로 부르지 않는다.

Bootstrap마다 50,000점의 예측을 다시 계산하지 않는다. 각 고정 cell/learner/route/rule 묶음에서 $G_{k\ell}\triangleq\mathbb P_T[f_kf_\ell]$, $v_k\triangleq\mathbb P_T[f_k\tau]$, $t\triangleq\mathbb P_T\tau^2$를 한 번 계산한다. 재표집 빈도에서 얻은 가중치 $a_k$의 합을 1로 두면 평균 위험은 $\sum_k a_k(G_{kk}-2v_k+t)$, 제곱편향은 $a^\top Ga-2a^\top v+t$, 경험분산 성분은 $\sum_k a_kG_{kk}-a^\top Ga$다. 기존 집계 함수 안에서 이 연산을 재사용하여 추가 모형 적합과 대형 grid 반복 연산을 피한다.

추가 적합 없이 고정 $(0.5,0.5)$ 후보를 진단한다. 이 후보와 RCT-only의 selection 점수 차이에서 모집단 위험 차이를 빼면 공통 noise floor가 상쇄된다. 정확 비율에서 그 잔차의 조건부 기댓값은 0이고, 추정·wrong-unit에서는 원고의 $\omega\mathbb E_O[(r-r_0)(\widehat g-\widehat\zeta)^2]$가 남는다. 실제 계산은 유한 50,000점 truth grid의 근사 위험을 사용하므로 적분 오차가 더해진다. 따라서 고정 grid에 조건부인 실제 잔차의 기댓값이 정확히 0이라고 주장하지 않으며 이 값은 정리 검증이 아닌 진단이다. 임의 절대 PASS cutoff를 두지 않는다. 같은 OBS selection 배열에서 ratio 평균·제곱평균·ESS·최댓값과 nuisance classifier 수렴을 기록한다. 표본 최댓값을 정리의 전역 상한으로 쓰지 않는다.

1. **그림 1, 위험 비교:** DRF/RF를 행, shared-unit / shifted-exact / shifted-classifier / shifted-wrong-unit을 열로 둔다. 네 규칙의 평균 위험비와 paired 구간, 기준선 1을 표시한다. Learner별 RCT baseline을 혼용하지 않는다.
2. **그림 2, 위험 분해:** 같은 패널·규칙에서 $\widehat B^2$와 $(K-1)\widehat V/K$를 쌓고 총높이 $\overline R$를 표시한다. 패널의 동일 RCT 평균 위험으로 두 성분을 함께 정규화한다.
3. **결과표 1:** cell / learner / ratio / rule / $K$ / mean risk / paired difference·구간 / variance / squared bias / failures를 기록한다. 절대값을 보존하고 regret·선택계수·고정후보 drift·ratio 진단은 보조 열로 둔다.

정확 비율이 복원하는 것은 기대 손실의 목표다. 무이동과 동일한 유한표본 위험, ratio 추정의 무비용, 모든 반복의 개선을 뜻하지 않는다. 수치실험은 이 SCM·예산·기저의 유한표본 결과를 보여준다. 기대 항등식의 정확성은 원고와 대수 검산에 근거한다.

## 8. 필요한 검증과 종료 기준

| 검사 | PASS 기준 |
| --- | --- |
| 의사결과 | 고정 배열에서 AIPW 직접식과 `build_scores` 일치 |
| 두 적합식 | 원래 목적함수 직접 미분과 정규방정식 일치; 상대 잔차 $10^{-10}$ 이하 |
| RF 가중치 | $e=0.5$에서 DRF/RF 일치, $e=0.35$에서 $\chi$ 누락 mutation 검출 |
| 선택식 | OBS 보정 제거·별도 omega·RF selection chi 삽입 mutation 검출 |
| 역할·누출 | 역할 disjoint; truth tau·reporting 값을 바꿔도 후보와 선택 불변 |
| 목표·ratio | RCT 설정 동일; Gaussian log-density 차이와 analytic ratio 일치; 이동 0이면 ratio 1 |
| 공정한 후보 | 같은 격자·ridge·자료·tie, omega0 경로 간 동일, true risk 선택은 oracle 열만 |
| 함수 일관성 | 적합·선택·truth에서 같은 함수, 숨은 clipping·penalty 없음 |
| 분해·집계 | 동일 grid 분해 상대오차 $10^{-10}$ 이하, paired bootstrap 단위·집계 mutation 검출 |
| 회계·재현 | 예정=성공+실패, seed 교체 없음, 설정·hash·원시예측에서 표·그림 재생성 |

구현 완료는 위 검사와 smoke 통과다. 실험 완료는 main 회계와 그림 2개·표 1개가 재현되는 상태다. 수치 실패와 불리한 결과도 보고한다. 식·누출·목표법칙 오류이면 필요한 부분만 수정하고 오염된 실행 묶음을 새 버전으로 다시 수행한다. 좋은 성능을 PASS 조건으로 두지 않는다.

기존 저장 방식에 새 run ID를 사용하여 설정·seed·source/code hash·raw 후보 점수·선택 결과·truth 예측·summary·failures·runtime·그림·표를 보존한다. 과거 결과는 덮어쓰지 않는다. 새 promotion framework는 필요 없다. 실제 재현 명령과 run 위치는 후속 구현 완료 시 기록한다.

## 9. 현재 제외한 항목

STAR 실제 결과 재구성, NSW, WHI, neural·별도 linear 기저, Riesz·balancing·hybrid 비교, 표본크기·여러 SCM·교란 축, 네 사분면, oracle 계수 회수·적응 Monte Carlo, 유계 정리 전용 DGP·반경 실험, ATE 재실행과 manuscript 수정은 보류한다. 현재 두 질문을 답하는 데 필요하지 않다. 추가하려면 어느 미해결 원고 주장에 필요한지와 기존 예측으로 답할 수 없는 이유부터 적는다.

---

<details>
<summary>비활성 역사: 2026-09-12 구 실행계획. 아래 요구·명령·승인 순서는 현재 적용하지 않는다.</summary>

# 구 DataFusionPPI CATE 확인 실험 실행 PRD

작성일: 2026-09-12

상태: 실행 전 operational contract

주대상: CATE confirmatory execution

## 0. Output Contract

### 0.1 Deliverable

이 문서는 DataFusionPPI의 CATE 실험을 다음 순서로 실행하기 위한 제품 요구사항 문서다.

1. 현재 M0/M1 인프라의 재현 확인
2. confirmatory validity blocker 수정과 재검증
3. M2 end-to-end smoke
4. M3 미니실험과 M4 의사결정용 확대
5. M5 synthetic·shift 확인 실험
6. repaired STAR 실험
7. 결과 집계, 그림, 표, manuscript handoff

이 문서는 통계 공식을 새로 정하지 않는다. 이미 정한 공식을 실제 명령, 입력, 반복
단위, 산출물, 중단 조건에 연결한다.

### 0.2 Consumer

용한은 이 문서로 다음 실행을 승인하거나 중단한다. 후속 구현자는 각 Stage의
`Purpose`, `Prerequisites`, `Command`, `Schedule`, `Outputs`, `PASS/FAIL`,
`Stop/Resume/Rollback`을 실행 계약으로 사용한다.

### 0.3 Completion criteria

전체 CATE 확인 실험이 완료되었다고 말하려면 다음 조건을 모두 만족해야 한다.

- 현 code fingerprint에서 M0/M1 authority가 존재한다.
- P0 blocker가 자동 테스트와 M2 smoke를 통과한다.
- M3 뒤, M4 전에 scientific decision record가 동결된다.
- `GO_CONFIRMATORY`가 기록된 경우에만 M5를 실행한다.
- 모든 예정 seed가 성공 결과 또는 명시적 실패 결과를 갖는다.
- synthetic 결과와 real-data reporting score를 구분한다.
- raw, summary, plot, table, provenance가 서로 hash로 연결된다.
- negative result와 tail failure를 삭제하지 않는다.
- manuscript에 들어갈 문장은 claim-to-evidence map을 통과한다.

### 0.4 Boundaries

이 문서 작성은 실험 실행이나 코드 변경 승인이 아니다. 실제 구현, 실행, manuscript
수정, mirror sync, stage, commit, push는 각각 승인된 범위에서만 한다. ATE는 이미 완료된
별도 evidence로 취급하며 이 PRD에서 재실행하지 않는다.

## 1. Authority hierarchy

충돌이 생기면 다음 순서로 판단한다.

1. 통계 식과 theorem 조건: `2026-09-09-implementation-handoff-JA.md` version 3
2. 논문의 claim과 estimator 정의: `manuscript/main/4.tex`
3. CATE 과학 설계와 감사 결과:
   [`2026-09-12-cate-experiment-audit-and-confirmatory-plan.md`](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md)
4. 실행 순서, 파일 계약, 복구 규칙: 현재 PRD
5. 한 stage의 실제 schedule과 seed: 해당 stage의 frozen `protocol.json`과
   `seed_ledger.json`

현재 PRD가 handoff의 식이나 theorem 적용 범위를 바꾸면 handoff가 우선한다. 기존 감사
문서는 현재 evidence와 새 DGP의 과학적 이유를 설명한다. 현재 PRD는 그 내용을 반복하지
않고 실행 가능한 단위로 바꾼다.

근거:

- 공식 CATE target과 risk: [audit plan §2.1](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#21-cate-estimand)
- manuscript-faithful loss와 oracle: [audit plan §6](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#6-manuscript-faithful-method-target)
- P0 요구: [audit plan §8](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#8-implementation-prd)
- M0부터 M5 설계: [audit plan §9](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#9-agile-execution-m0-to-m5)

## 2. 처음 쓰는 용어와 평가량

### 2.1 CATE와 true risk

$X$는 처치 전에 측정한 공변량이다. $Y(1)$과 $Y(0)$은 각각 처치와 대조 상태의
potential outcome이다. RCT target population의 조건부 평균 처치효과를

$$
\tau(x)=\mathbb E_R\{Y(1)-Y(0)\mid X=x\}
$$

로 쓴다. 이것이 CATE다. 추정 함수 $\widehat\zeta$의 integrated squared risk는

$$
\mathcal R(\widehat\zeta)
=\mathbb E_R[\{\widehat\zeta(X)-\tau(X)\}^2]
$$

다. Synthetic DGP에서는 $\tau$를 아니까 independent truth sample로 계산한다. STAR,
NSW, WHI에서는 pointwise $\tau$를 모르므로 이 값을 보고하지 않는다
([audit plan:116](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#21-cate-estimand)).

### 2.2 네 data role

한 source sample이 다음 네 역할을 겸하지 않게 한다.

- **nuisance:** outcome regression, propensity, representation을 학습한다.
- **tuning:** 각 candidate coefficient를 적합한다.
- **selection:** candidate 중 하나를 선택한다.
- **reporting:** 선택이 끝난 뒤 성능 score를 계산한다.

Selection이 reporting outcome을 보면 선택 후 score가 낙관적으로 바뀔 수 있다. 그래서
두 역할은 독립이어야 한다. 현재 exploratory engine은 세 역할을 사용하므로 P0에서
고쳐야 한다
([audit plan:206](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#31-공통-engine),
[audit plan:312](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#42-honest-selection-and-reporting)).

### 2.3 두 channel과 coefficient

$\lambda$는 RCT와 OBS outcome regression을 섞어 RCT-valid score의 noise를 줄이는
channel이다. $\omega$는 OBS에서 학습한 CATE prediction과 큰 OBS covariate sample을
이용해 fitting variance를 줄이는 channel이다. $\rho$는 linear basis coefficient에 주는
ridge penalty다. Candidate는 $(\lambda,\omega,\rho)$의 한 조합이다.

### 2.4 Transport ratio

$P_R^X$와 $P_O^X$는 RCT와 OBS의 $X$ 분포다. 두 분포가 다르면

$$
r_0(x)=\frac{dP_R^X}{dP_O^X}(x)
$$

로 OBS average를 RCT target으로 옮긴다. Gaussian shift에서는 exact ratio가 있어도
unbounded다. STAR primary source law는 $P_O^X=P_R^X$를 구성하여 $r_0=1$로 만든다.
이 equality는 causal identification을 증명하지 않는다
([audit plan §2.3](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#23-transport-ratio)).

### 2.5 서로 바꾸어 쓸 수 없는 수치

- **True risk:** synthetic에서만 계산하는 $\mathcal R(\widehat\zeta)$.
- **Reporting-score difference:** real data의 independent RCT reporting sample에서 계산한
  selected minus RCT-only score. True CATE risk가 아니다.
- **Prediction variance:** 같은 target point에서 replication에 따라 prediction이 변하는
  정도를 적분한 값.
- **Squared bias:** replication 평균 prediction과 $\tau$의 차이를 제곱해 적분한 값.
- **Selection regret:** selected risk minus 같은 grid의 true-risk oracle risk.
- **Risk ratio:** selected risk divided by matched RCT-only risk. 분모가 작으면 불안정하다.
- **Failure rate:** 예정 unit 중 numerical, protocol, scientific failure의 비율.

Independent replication 수는 서로 독립으로 새 data를 생성한 횟수다. 한 replication에서
learner, route, candidate, metric 행이 많이 생겨도 반복 수는 늘지 않는다.

## 3. Current-state distinction

### 3.1 공식 v3 G0

공식 v3 G0는 ATE와 CATE가 섞인 8개 performance cell에서 replication ID 0을 한 번씩
실행했다. 독립 반복은 8회이고 372는 long-format metric row 수다. 결정은
`INCONCLUSIVE`다
([implementation roadmap:61](./implementation-experiment-roadmap.md#31-완료된-g0)).

### 3.2 CATE confirm M0/M1

`materials/confirm`의 현재 promoted M0/M1은 모두 historical/pre-fix다. 현재 code
fingerprint의 authority는 아직 없다
([confirm README](./confirm/README.md#historical-pre-fix-stages)).

현재 `code/confirm/m1_protocol.py`가 정한 M1 계약은 다음과 같다.

- scheduled units: 26
- expected checks: 53
- required promoted outputs: 8
- seed key: protocol ID, stage, base seed, unit

근거: [m1_protocol.py](../code/confirm/m1_protocol.py).

### 3.3 Existing CATE results

CATE-1, CATE-2, CATE-3 CSV와 그림은 exploratory evidence다. 그 결과는 DGP를 개선하고
failure mode를 찾는 데 사용하지만 confirmatory 결과와 합치지 않는다. 현재 CATE-1은
five-family exploratory grid이고, CATE-2 STAR는 v3 exact source law가 아니며, CATE-3는
estimated-ratio extension이다
([audit plan §3](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#3-현재-implementation-inventory)).

## 4. 공통 artifact와 provenance 계약

### 4.1 Planned artifact tree

각 정상 실행은 다음 형태의 새 디렉터리 하나를 만든다.

```text
materials/confirm/
  <STAGE>-<code-fingerprint-prefix>/
    FINGERPRINT
    code_manifest.json
    protocol.json
    seed_ledger.json
    schedule.json
    environment.json
    checks.json
    failures.json
    runtime.json
    RECORD.md
    raw/
      replications-<shard>.csv
    summary/
      cell_summary.csv
      mechanism_summary.csv
      tail_summary.csv
      ratio_summary.csv
    figures/
    tables/
```

M0와 M1은 이미 정한 더 작은 exact output set을 그대로 사용한다. M2 이후 stage만 위
공통 tree를 사용한다.

### 4.2 Required raw row keys

모든 synthetic raw row에는 최소한 다음 열이 있어야 한다.

```text
protocol_id, code_fingerprint, stage, study, family, cell_id,
replication_id, seed_key, seed32, role_fingerprints,
learner, sieve, ratio_route, candidate_id, lambda, omega, rho,
metric, value, status, failure_class, runtime_seconds
```

STAR row에는 `outcome`, `alpha`, `source_law_fingerprint`, `reporting_score`를 더한다.
Synthetic 전용 `true_risk`는 STAR에서 `N/A`여야 한다.

### 4.3 Fingerprint와 no-overwrite

- default stage name은 full code fingerprint의 prefix를 포함한다.
- explicit stage name은 preflight에서 collision 상태를 보여준다.
- existing promoted stage를 다른 fingerprint로 교체하지 않는다.
- 같은 fingerprint 재실행도 명시적 replace 승인 없이는 거부한다.
- raw shard는 immutable이다.
- summary와 plot은 raw hash를 input manifest에 기록한다.
- code fingerprint는 calculation 직전과 promotion 직전에 다시 확인한다.

현재 atomic promotion 구현 근거는
[provenance.py](../code/confirm/provenance.py)다.

### 4.4 Seed, failure, resume

- seed는 실행 전에 ledger에 존재해야 한다.
- 실패한 seed를 새 seed로 바꾸지 않는다.
- recoverable model failure는 raw failure row로 남긴다.
- hard provenance 또는 schema failure는 promotion을 중단한다.
- resume은 ledger에서 완료되지 않은 unit만 새 scratch에서 계산한다.
- resume 결과는 기존 immutable shard와 hash가 맞을 때만 함께 promotion한다.
- scheduled count는 success count와 explicit failure count의 합과 같아야 한다.

### 4.5 Rollback

Rollback은 파일 삭제가 아니다. 새 stage promotion이 실패하면 이전 promoted authority를
계속 사용한다. 현재 Stage instance가 만든 scratch만 지운다. 다른 process의 scratch와
historical debris는 scan하거나 삭제하지 않는다.

## 5. Current CLI와 planned CLI

### 5.1 현재 존재하는 명령

```bash
cd /Users/yonghanjung/paios/research/papers/DataFusionPPI/code
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest test_v3_gate.ConfirmInfrastructureTests
PYTHONDONTWRITEBYTECODE=1 python3 -m confirm.m0_freeze --preflight
PYTHONDONTWRITEBYTECODE=1 python3 -m confirm.m1_fixtures --preflight
python3 -m confirm.m0_freeze
python3 -m confirm.m1_fixtures
```

마지막 두 명령은 실제 stage를 쓰므로 별도 실행 승인이 필요하다.

### 5.2 앞으로 구현할 명령

다음 명령은 아직 존재하지 않는다. 문서의 예시는 planned CLI contract다.

```bash
python3 -m confirm.cate_run --stage m2 --preflight
python3 -m confirm.cate_run --stage m2
python3 -m confirm.cate_run --stage m3 --preflight
python3 -m confirm.cate_run --stage m3
python3 -m confirm.cate_run --stage m4 --decision-record <path>
python3 -m confirm.cate_run --stage m5 --decision-record <path>
python3 -m confirm.cate_run --stage star --decision-record <path>
python3 -m confirm.cate_analyze --input-stage <path> --preflight
python3 -m confirm.cate_analyze --input-stage <path>
python3 -m confirm.cate_plot --analysis-stage <path> --preflight
python3 -m confirm.cate_plot --analysis-stage <path>
```

계획 구현 파일은 다음 여섯 개가 최소 단위다. 생성 또는 수정은 별도 승인을 받는다.

- `code/confirm/cate_protocol.py`
- `code/confirm/cate_engine.py`
- `code/confirm/cate_run.py`
- `code/confirm/cate_analyze.py`
- `code/confirm/cate_plot.py`
- `code/test_cate_confirm.py`

기존 exploratory driver를 조용히 confirmatory driver로 바꾸지 않는다.

## 6. Stage 0: Read-only preflight

### Purpose

실행 전에 source, target, schedule, collision을 확인한다. 과학 계산과 결과 쓰기는 하지
않는다.

### Prerequisites

- canonical project가 접근 가능함
- required input 문서와 code가 존재함
- 다른 process가 target stage를 promotion 중이지 않음

### Exact command

현재 명령:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m confirm.m0_freeze --preflight
PYTHONDONTWRITEBYTECODE=1 python3 -m confirm.m1_fixtures --preflight
```

M2 이후에는 같은 contract의 planned `confirm.cate_run --stage <stage> --preflight`를 쓴다.

### Scheduled units and computation

0 scientific units, 0 model fits, 0 result writes.

### Required output

stdout의 machine-readable JSON만 허용한다.

- `calculation_performed=false`
- `write_performed=false`
- resolved target
- code fingerprint
- protocol and schedule hash
- collision state
- required input/output 목록
- qualification status

### PASS/FAIL

PASS는 모든 input이 존재하고, schedule이 exact하며, target이 qualified되고, collision이
설명 가능한 경우다. Missing input, schedule mismatch, unqualified default, 설명되지 않은
collision은 FAIL이다.

### Stop/Resume/Rollback

FAIL이면 정상 명령을 실행하지 않는다. Preflight는 stage를 만들지 않으므로 rollback이
필요 없다.

## 7. Stage 1: Current-code M0/M1 reproduction

### Purpose

P0 수정 전 코드 $F_0$의 evidence와 deterministic mechanism 상태를 재현 가능한 authority로
남긴다. $F_0$는 이 단계에서 계산한 full code fingerprint다.

### Prerequisites

- Stage 0 PASS
- infrastructure tests PASS
- explicit write approval
- historical M0/M1을 current라고 부르지 않음

### Exact command

현재 구현된 명령:

```bash
python3 -m confirm.m0_freeze
python3 -m confirm.m1_fixtures
```

### Scheduled units and computation

- M0: artifact, code, document freeze와 promotion self-test
- M1: 4 base + 16 mutations + 4 constant-propensity RF + sampler + bootstrap
- 합계: 26 deterministic/stochastic fixture units
- 성능 replication: 0

M1 schedule 근거는 [m1_protocol.py](../code/confirm/m1_protocol.py)다.

### Required outputs

M0 exact 8:

```text
FINGERPRINT, artifact_manifest.json, code_manifest.json,
document_manifest.json, seed_ledger.json, environment.json,
label_map.json, promotion_self_test.txt
```

M1 exact 8:

```text
FINGERPRINT, code_manifest.json, protocol.json, seed_ledger.json,
environment.json, checks.json, moments.json, RECORD.md
```

### PASS/FAIL

- 26-unit ledger exact, seed collision 0
- 53/53 check PASS
- production과 독립 moment formula 차이 $\le10^{-10}$
- required output의 누락과 extra file 0
- code fingerprint가 시작과 promotion 시점에 같음

하나라도 어기면 FAIL이다.

### Stop/Resume/Rollback

M1 FAIL이면 P0 구현과 원인 분석은 가능하지만 M2 실행은 금지한다. Failed scratch를
promoted stage처럼 쓰지 않는다. $F_0$ stage는 P0 뒤 결과와 합치지 않는다.

## 8. Stage 2: P0 blocker implementation and post-fix gate

### Purpose

성능을 보기 전에 confirmatory validity를 확보한다. P0 수정 뒤 code fingerprint를 $F_1$로
정의한다.

### Prerequisites

- Stage 1의 재현 결과 또는 Stage 1을 생략한 이유를 decision record에 기록
- P0 구현 범위와 파일 목록 승인
- 기존 exploratory raw 결과를 수정하지 않음

### P0 blockers

| blocker | 구현 요구 | 주요 tests |
| --- | --- | --- |
| four-role engine | source별 nuisance, tuning, selection, reporting 독립 | T04-T07 |
| theorem paths | bounded Theorem 6과 unclipped $\rho=0$ Theorem 7 분리 | T12-T14 |
| adaptive oracle | 20k, 40k, 80k, 최대 100k와 moment별 MCSE | T15-T18 |
| STAR source law | exact $X$, OOF rank, exact-cell propensity, frozen law | T08-T11, T33 |
| failure and tail ledger | every seed accounting, no redraw, tail classification | T27-T29 |
| endpoints | paired risk, decomposition, calibration, ratio diagnostics | T19-T26, T30-T32 |

자세한 과학 조건은 [audit plan §8.1](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#81-p0-confirmatory-validity)을 따른다.

### Planned command

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v test_cate_confirm.py
python3 -m confirm.m0_freeze --preflight
python3 -m confirm.m1_fixtures --preflight
python3 -m confirm.m0_freeze
python3 -m confirm.m1_fixtures
```

첫 줄과 post-fix M0/M1 normal run은 계획 명령이다. 실제 파일 이름과 test class는 구현
승인 때 확정한다.

### Scheduled units and computation

- fixed-array unit tests와 small deterministic fixtures만 실행
- performance replication 0
- post-fix M1은 같은 26-unit contract를 새 $F_1$ stage에서 재실행

### Required outputs

- test report with T01-T35 status
- P0 implementation manifest
- $F_0$ to $F_1$ correspondence map
- post-fix M0/M1 exact outputs
- unresolved test와 explicit waiver 0

### PASS/FAIL

Applicable hard tests 전부 PASS해야 한다. Reporting outcome perturbation은 candidate와
coefficient를 바꾸지 않아야 한다. Theorem 7 table에는 neural, clipped spline,
$\rho>0$ row가 없어야 한다. STAR source-law와 exact-cell propensity 오차는 각각
$10^{-12}$ 이하여야 한다.

### Stop/Resume/Rollback

P0 hard failure는 `REDESIGN`이다. $F_0$ M0/M1을 $F_1$ 결과 prefix로 자동 재사용하지
않는다. 수정 일부를 되돌릴 때도 historical stage를 삭제하지 않고 새 fingerprint를 쓴다.

## 9. Stage 3: M2 one-rep end-to-end smoke

### Purpose

작은 비용으로 data generation부터 immutable raw row, validation, atomic promotion까지
전체 pipeline이 한 번 연결되는지 확인한다.

### Prerequisites

- Stage 2 PASS
- M2 schedule, seed ledger, result schema 동결
- planned runner와 validator unit tests PASS
- M2 stage preflight PASS

### Planned command

```bash
python3 -m confirm.cate_run --stage m2 --preflight
python3 -m confirm.cate_run --stage m2
python3 -m confirm.cate_run --stage m2 --validate-only <M2-stage-path>
```

### Exact proposed schedule

| unit | source unit 수 | 같은 source unit에서 실행할 path |
| --- | ---: | --- |
| L0O0 spline DRF | 1 | RCT, $\lambda$, $\omega$, joint, grid oracle |
| L1O0 spline DRF | 1 | 같은 estimator set |
| L0O1 spline DRF | 1 | 같은 estimator set |
| L1O1 spline DRF | 1 | 같은 estimator set |
| Gaussian strong shift | 1 | oracle, classifier, BAL-X, BAL-X+$g$, unit |
| STAR math, $\alpha=0.8,n_R=200$ | 1 | frozen source law, select, independent report |
| 합계 | 6 | 6 independent stochastic source units |

다음 검사는 독립 source unit를 추가하지 않는다.

- constant-propensity RF cross term
- bounded Theorem 6와 `bound_vacuity_ratio`
- unclipped $\rho=0$ Theorem 7
- neural empirical path
- ratio balance·normalization·ESS
- STAR role and source-law assertions

이 6-unit proposal은 실제 M2 구현 전에 `cate_protocol.py`에 동결하고 통계 리뷰를 받는다.
Route나 candidate row 수를 6에 곱해 반복 수라고 부르지 않는다.

### Required outputs

- exact six-unit `schedule.json`
- six successful or explicit-failure source records
- long-format raw rows
- role and source fingerprints
- theorem eligibility flags
- ratio diagnostics
- failure and runtime ledgers
- `checks.json`과 deterministic `RECORD.md`

### PASS/FAIL

- 모든 P0 hard assertion PASS
- scheduled = success + explicit failure
- synthetic true-risk fields finite 또는 명시적 failure
- STAR true-risk field `N/A`
- five ratio routes exact
- result schema와 fingerprint exact

Joint가 RCT-only보다 나쁜 방향은 engineering FAIL이 아니다. Formula, role, source law,
propensity, theorem label, provenance failure는 `REDESIGN`이다.

### Stop/Resume/Rollback

Hard failure이면 M3를 잠근다. Interrupted M2는 completed unit hash를 검증한 뒤 missing
unit만 resume한다. Partial scratch는 evidence가 아니다.

## 10. Stage 4: M3 ten-rep mini experiment

### Purpose

Finite output, 방향, 효과 크기, tail mechanism, runtime을 측정한다. Superiority를
확정하지 않는다.

### Prerequisites

- M2 PASS
- paired analysis와 bootstrap fixture PASS
- M3 schedule과 candidate grid 동결
- outcome을 보기 전에 aggregation code hash 동결

### Planned command

```bash
python3 -m confirm.cate_run --stage m3 --preflight
python3 -m confirm.cate_run --stage m3
python3 -m confirm.cate_analyze --input-stage <M3-stage-path>
```

### Exact schedule

| panel | cells | reps/cell | independent reps |
| --- | ---: | ---: | ---: |
| 4 quadrants × 2 confounding shapes × 3 $n_R$ | 24 | 10 | 240 |
| O1/L1O1 × 2 shifts × 3 $n_R$ | 12 | 10 | 120 |
| 6 negative/failure controls | 6 | 10 | 60 |
| 합계 | 42 |  | 420 |

근거: [audit plan:819](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#m3-ten-rep-pilot).

### Required outputs

- 420 source-unit accounting
- cell·learner·sieve·route runtime
- primary paired differences의 pilot variance
- risk and ratio tail quantiles
- failure categories and seeds
- score-risk calibration pilot
- M4/M5 runtime forecast input

### PASS/FAIL

M3의 PASS는 우월성이 아니라 pipeline과 planning information의 완전성이다. 다음이 있어야
한다.

- scheduled accounting 100%
- finite/failed status 100%
- paired keys exact
- runtime and peak memory recorded
- catastrophic event가 발생하면 seed와 raw row 보존

### Stop/Resume/Rollback

Hard validity failure이면 `REDESIGN`이다. Scientific direction이 약하면 결과를 숨기지 않고
M4 decision input으로 쓴다. 현재 sequence에서 M3는 항상 exploratory design pilot이다.
DGP, grid, threshold를 바꾸지 않더라도 M5 primary inference에 재사용하지 않는다.

## 11. Stage 5: Decision record and M4 escalation

### Purpose

아직 불확실한 cell만 늘리고, M5를 열지 여부를 한 개의 machine-readable decision으로
기록한다.

### Prerequisites

- M3 raw와 summary freeze
- M3 runtime report
- scientific threshold를 정할 책임자와 시점 확정

### Decision record에서 동결할 미정 항목

현재 다음 값은 확정하지 않는다.

- minimum material effect
- noninferiority margin
- interval confidence level과 bootstrap resample count
- family aggregation과 multiplicity adjustment
- catastrophic-risk cutoff
- acceptable scientific negative-transfer rate

이 값은 M3가 끝난 뒤, M4 outcome을 보기 전에 근거와 함께 동결한다.
기존 5% 기준은 protocol choice이지 통계 법칙이 아니다
([audit plan:835](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#m3-ten-rep-pilot)).

### Planned command

```bash
python3 -m confirm.cate_run --stage m4 --decision-record <frozen-json> --preflight
python3 -m confirm.cate_run --stage m4 --decision-record <frozen-json>
```

### Schedule and cap

- M3의 42 cells 중 relevant inconclusive cell만 증가
- cell별 누적 cap 50
- 전 cell을 채울 경우 누적 최대 2,100
- M3 420회를 포함하므로 추가 최대 1,680
- 자동 full-cap 실행 금지

### Required outputs

- immutable decision JSON
- increment schedule and rationale
- cumulative accounting
- updated interval and tail summaries
- terminal decision exactly one

허용 terminal token:

- `REDESIGN`
- `REFRAME`
- `INCONCLUSIVE`
- `GO_CONFIRMATORY`

### PASS/FAIL

- hard validity failure: `REDESIGN`
- valid하지만 channel opportunity 없음: `REFRAME`
- valid하고 cap 아래 uncertainty가 큼: `INCONCLUSIVE`
- prespecified scientific and safety criteria 충족: `GO_CONFIRMATORY`

### Stop/Resume/Rollback

`GO_CONFIRMATORY` 외 상태는 M5를 잠근다. M4 중 protocol을 바꾸면 새 version을 만들고
기존 결과를 섞지 않는다. M3와 M4는 모두 exploratory design data이며 M5 primary
estimate나 confidence interval에 합치지 않는다.

## 12. Stage 6: M5 full synthetic and shift simulation

### Purpose

FusionPPI가 도움되는 조건, 중립인 조건, 위험해지는 조건을 같은 frozen protocol에서
확인한다.

### Prerequisites

- M4 decision이 정확히 `GO_CONFIRMATORY`
- final protocol, grid, DGP, endpoints, multiplicity, tail threshold freeze
- M3/M4와 겹치지 않는 confirmatory seed namespace freeze
- runtime and storage budget 승인
- M5 preflight PASS

### Planned command

```bash
python3 -m confirm.cate_run --stage m5 --decision-record <frozen-json> --preflight
python3 -m confirm.cate_run --stage m5 --decision-record <frozen-json>
python3 -m confirm.cate_run --stage m5 --validate-only <M5-stage-path>
```

### Exact schedule

| panel | cells | reps/cell | independent reps |
| --- | ---: | ---: | ---: |
| 4 quadrants × 2 confounding shapes × 3 $n_R$ | 24 | 200 | 4,800 |
| O1/L1O1 × 2 shifts × 3 $n_R$ | 12 | 200 | 2,400 |
| 6 negative/failure controls | 6 | 100 | 600 |
| 합계 | 42 |  | 7,800 |

근거: [audit plan:853](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#m5-confirmatory-simulation).

7,800은 unique DGP replications다. Model fit 수는 learner, sieve, route, candidate를 반영해
별도 runtime manifest에서 계산한다.

### Required outputs

- 7,800-unit raw or failure ledger
- M3/M4와 disjoint한 seed namespace 검증
- cell summaries with paired uncertainty
- mechanism moments and coefficients
- variance/bias decomposition
- calibration and regret
- ratio and overlap diagnostics
- tail and runtime summaries
- exact raw and summary manifests

### PASS/FAIL

Execution PASS는 결과 방향과 무관하게 schema, accounting, provenance, prespecified analysis가
완전한 상태다. Scientific claim은 frozen decision criteria에 따라 별도로 PASS, REFRAME,
negative result로 기록한다.

### Stop/Resume/Rollback

Hard failure는 새 calculation을 중단한다. 이미 완료된 immutable shard는 보존한다. Code
또는 protocol drift가 생기면 새 fingerprint stage로 다시 시작한다. 실패 seed를 교체하지
않는다.

## 13. Stage 7: Repaired STAR full run

### Purpose

실제 STAR covariate와 outcome을 사용하는 source-marginal stress test를 수행한다. True CATE
recovery 실험이 아니다.

### Prerequisites

- P0 STAR exact source-law tests PASS
- data source, license, SHA-256, schema freeze
- eligible exact-$X$ cells의 $e_x\in[0.15,0.85]$
- OOF residual rank와 source weights freeze
- learner features에서 rank, residual, weight 제외

### Planned command

```bash
python3 -m confirm.cate_run --stage star --decision-record <frozen-json> --preflight
python3 -m confirm.cate_run --stage star --decision-record <frozen-json>
```

### Exact schedule

| axis | values |
| --- | --- |
| outcome | mathematics, reading |
| RCT role size $n_R$ | 200, 400 |
| selection strength $\alpha$ | 0.4, 0.8 |
| replications | 100 per cell |

총 2×2×2×100=800 independent source-law replications다.

Role sizes:

- RCT nuisance and tuning: 각각 $n_R$
- OBS nuisance and tuning: 각각 5,000
- selection and reporting: source별 각각 1,000

근거: [audit plan §11.1](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#111-repaired-star-primary-stress-test).

### Required outputs

- frozen support, fold, probability, propensity fingerprints
- source-level exact marginal and propensity checks
- independent reporting-score differences
- eligible-cohort descriptive ATE benchmark
- prediction and coefficient stability
- failure and runtime ledger

### PASS/FAIL

- $\max_x|P_O^X(x)-P_R^X(x)|\le10^{-12}$
- $\max_x|e(x)-P_R(A=1\mid X=x)|\le10^{-12}$
- role RNG independence
- reporting perturbation이 selection을 바꾸지 않음
- true-risk field `N/A`

Score 방향이 불리한 것은 code FAIL이 아니다. Source law, role, propensity, leakage,
provenance failure는 hard FAIL이다.

### Stop/Resume/Rollback

Source-law freeze 실패 시 STAR 전체를 중단한다. Realized sample의 $X$ 비율이 우연히 다른
것은 source-law failure가 아니다. 동일 frozen law와 seed ledger에서만 resume한다.

## 14. Stage 8: Supporting real-world studies

### 14.1 NSW with CPS and PSID

NSW experimental treated/control을 RCT source로, CPS와 PSID comparison controls를 각각 OBS
source로 사용한다. 두 comparison은 합치지 않고 별도 panel로 보고한다. Primary
confirmatory real-data claim이 아니라 supporting benchmark다
([audit plan §11.2](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#112-nsw-with-cps-and-psid-supporting-benchmark-only)).

실행 전 다음을 동결한다.

- 1978 earnings estimand
- exact variable crosswalk
- experimental target law와 comparison law의 ratio target
- treated unit 중복을 막는 role partition
- source URL, version, hash, redistribution terms

평가는 independent experimental reporting score와 experimental ATE benchmark다. Pointwise
CATE truth나 oracle ratio를 만들지 않는다.

### 14.2 WHI feasibility-gated optional study

WHI Clinical Trial과 Observational Study는 자연스러운 RCT/OBS pair 후보지만 다음 gate를
먼저 통과해야 한다.

- participant-level access와 DUA
- CT/OS common-variable crosswalk
- treatment initiation와 follow-up alignment
- event count, overlap, censoring support
- 공개 가능한 artifact 범위

Binary 5-year risk difference와 restricted mean survival time 중 하나를 outcome inspection
전에 고정한다. Survival estimand를 고르면 censoring-aware loss를 먼저 구현한다. Feasibility
FAIL이면 WHI를 실행하지 않는다
([audit plan §11.3](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#113-whi-ct-plus-os-strongest-natural-pair-feasibility)).

## 15. Stage 9: Aggregation and uncertainty

### Purpose

Frozen raw 결과에서 통계량을 한 번 정의된 방식으로 계산한다. Aggregation code는 raw를
수정하지 않는다.

### Primary synthetic estimand

Replication $k$에서 Fusion과 matched RCT-only의 true risk 차이를

$$
d_k=\mathcal R_k(\widehat\zeta_{\mathrm{Fusion}})
-\mathcal R_k(\widehat\zeta_{\mathrm{RCT}})
$$

로 둔다. Cell의 primary point estimate는 $d_k$의 평균이다. Uncertainty unit은 paired
replication이다. Family pooled interval은 family-stratified paired bootstrap을 쓴다
([audit plan §10.1](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#101-primary-synthetic-endpoint)).

### Required secondary summaries

| quantity | aggregation | uncertainty unit | caveat |
| --- | --- | --- | --- |
| absolute risk | cell별 replication mean | replication | scale를 함께 보고 |
| risk ratio | per-rep ratio의 median과 tails | paired replication | 작은 RCT risk에서 불안정 |
| prediction variance | fixed target anchors에 적분 | replication | anchor hash 필요 |
| squared bias | replication mean prediction 대 truth | replication and anchor | synthetic only |
| selection regret | selected minus grid oracle | paired replication | $<-10^{-12}$이면 hard fail |
| coefficient stability | boundary, fallback, distribution | replication | oracle type 분리 |
| calibration | score-risk slope, intercept, rank | candidate within replication | select/report 분리 |
| ratio diagnostics | route별 $E_Or$, $E_Or^2$, ESS, balance, drift | replication | normalization은 diagnostic |
| tail | 90, 95, 99%, maximum, failure rate | replication | threshold 사전 동결 |
| runtime | median, 90%, max, peak memory | fit and replication | row count로 나누지 않음 |
| real-data score | selected minus RCT reporting score | source replication | true risk 아님 |

Prediction variance와 squared bias는 fixed anchor에서

$$
\overline R_K=\frac{K-1}{K}\widehat V_K+\widehat B_K^2
$$

를 만족해야 한다. $K$는 독립 replication 수다. Decomposition gap은 deterministic fixture에서
$10^{-10}$ 이하여야 한다
([audit plan §8.2](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#82-p1-scientific-endpoints)).

### PASS/FAIL

- route-row가 아니라 unique family-cell-rep key를 count함
- paired comparison의 두 row가 모두 존재함
- missing과 failure를 0으로 바꾸지 않음
- synthetic와 STAR metric namespace가 분리됨
- bootstrap seed와 resample count가 frozen decision record와 일치함

## 16. Stage 10: Figures, tables, integrated handoff

### 16.1 Figure contract

| figure | content | mandatory safeguards |
| --- | --- | --- |
| C1 | $\lambda_p^\star$ 대 $\omega_p^\star$ opportunity map | oracle MCSE, quadrant label |
| C2 | paired mean risk difference | DRF/RF, spline/neural, zero line, paired CI |
| C3 | prediction variance and squared bias | matched RCT normalization, raw total annotation |
| C4 | coefficient recovery and regret | grid oracle와 moment oracle 구분 |
| C5 | five transport routes | oracle/classifier/BAL-X/BAL-X+$g$/unit 모두 표시 |
| C6 | tail and failures | numerical/protocol/scientific failure 분리 |
| C7 | real-data reporting score | STAR math/read, NSW CPS/PSID 분리, true-risk 표기 금지 |

세부 내용의 authority는
[audit plan §13](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#13-plot-prd)이다.

모든 caption은 metric 정의, aggregation estimand, 독립 replication count, uncertainty unit,
protocol status를 포함한다.

### 16.2 Table contract

| table | required content |
| --- | --- |
| C1 | design, cell, role sizes, rep IDs, code/protocol/source hashes |
| C2 | paired risk, interval, absolute risk, median ratio, independent $K$ |
| C3 | $A_p,B_p,C_p,D_p$, MCSE, oracle and selected coefficients, ablations |
| C4 | variance, squared bias, total risk, anchor size, decomposition gap |
| C5 | score-risk calibration, optimism, regret |
| C6 | five-route ratio and overlap diagnostics |
| C7 | failure classes, denominator, rate, seeds, mechanism |
| C8 | real-data estimand, construction, reporting score, ATE, license, claim |

세부 내용의 authority는
[audit plan §14](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#14-table-prd)이다.

### 16.3 ATE-CATE integrated minimum table

ATE는 기존 artifact manifest로 고정하고 재실행하지 않는다. 통합 표에는 최소한 다음 열을
넣는다.

```text
study, estimand, target population, estimator, comparator,
metric definition, point estimate, interval, independent replication count,
role contract, protocol status, code hash, result hash, allowed claim
```

다음은 같은 열에서 숫자만 나란히 놓더라도 동일 metric으로 해석하지 않는다.

- ATE bias, variance, MSE/RMSE, coverage
- CATE integrated true risk, prediction variance, squared bias
- STAR/NSW/WHI reporting-score difference

### 16.4 Manuscript handoff

Manuscript를 직접 수정하기 전에 다음 packet을 만든다.

- claim-to-evidence matrix
- figure/table inventory and SHA-256 manifest
- safe and forbidden claims
- negative result and failure ledger
- exact reproduction commands
- exploratory versus confirmatory artifact map
- unresolved theorem extensions

CATE-3 estimated-ratio 결과는 이론이 확장되지 않는 한 exploratory로 표시한다. Neural은
Theorem 7 evidence로 쓰지 않는다. STAR source-law equality는 causal identification으로
해석하지 않는다
([audit plan §12](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#12-safe-and-forbidden-claims)).

## 17. T01 to T35 acceptance mapping

아래 35개 requirement 중 33개는 hard condition이고, 하나는 report condition, 하나는
diagnostic condition이다. 35개 모두 machine-readable row를 가져야 하므로 90% 이상을
자동 validator가 확인한다.

| ID | stage | machine-observable PASS |
| --- | --- | --- |
| T01 | P0/M1 | pseudo-outcome identity gap $\le10^{-10}$ |
| T02 | P0/M1 | DRF solve와 direct objective gap $\le10^{-10}$ |
| T03 | P0/M1 | RF solve와 direct objective gap $\le10^{-10}$ |
| T04 | P0/M2 | four role IDs 또는 draw fingerprints 독립 |
| T05 | P0/M2 | reporting perturbation 뒤 candidate와 coefficient 불변 |
| T06 | P0/M2 | selection perturbation 뒤 reporting sample 불변 |
| T07 | all | namespace seed stable and distinct |
| T08 | P0/STAR | source-law marginal error $\le10^{-12}$ |
| T09 | P0/STAR | exact-cell propensity identity error $\le10^{-12}$ |
| T10 | P0/STAR | OOF prediction에서 own row 제외 |
| T11 | P0/STAR | learner $X$에 residual, rank, weight 없음 |
| T12 | P0/M2 | theorem table에 eligible path만 존재 |
| T13 | P0/M2 | bounded raw components가 prespecified bound 만족 |
| T14 | M2/M5 | radius와 vacuity ratio finite and reported |
| T15 | M1/P0 | production과 independent moments gap $\le10^{-10}$ |
| T16 | P0/M2 | adaptive oracle draw sequence와 stopping exact |
| T17 | M2/M5 | exact-ratio oracle risk가 candidate minimum |
| T18 | M2/M5 | selection regret $\ge-10^{-12}$ |
| T19 | P0/M3 | risk decomposition gap $\le10^{-10}$ on fixture |
| T20 | P0/M3 | score-risk identity gap $\le10^{-10}$ on discrete support |
| T21 | M3/M5 | candidate score와 risk가 key로 완전 pairing |
| T22 | P0/M2 | feasible BAL-X residual $\le10^{-6}$ |
| T23 | P0/M2 | BAL-X+$g$가 declared $g$ feature를 포함 |
| T24 | M2/M5 | five routes exact, omission 0 |
| T25 | P0/M2 | Gaussian analytic moments and ESS identity |
| T26 | all | Gaussian radius `N/A`, sample-max bound claim 0 |
| T27 | all runs | scheduled = success + explicit failure |
| T28 | all runs | duplicate or replacement seed 0 |
| T29 | P0/M3 | known exploratory tail seed 재현 또는 변경 설명 |
| T30 | P0/M3 | paired bootstrap fixture와 implementation 일치 |
| T31 | all summaries | unique family-cell-rep count 사용 |
| T32 | M2 onward | runtime, memory, code hash 존재 |
| T33 | STAR | true-risk `N/A`, reporting-score label 존재 |
| T34 | shift/report | `exploratory estimated-ratio extension` label 존재 |
| T35 | report | interval을 만들지 않은 결과에 coverage claim 0 |

원 acceptance matrix는
[audit plan §15](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#15-acceptance-test-matrix)에 있다.

## 18. Runtime and capacity budget

### 18.1 Before measured pilot timing

다음은 성능 예측이 아니라 운영 checkpoint다.

| stage | checkpoint |
| --- | ---: |
| Stage 0 preflight | 5분 |
| M0/M1 reproduction | 15분 |
| P0 unit and fixture gate | 15분 |
| M2 six-unit smoke | 60분 |

Checkpoint를 넘으면 중단하고 completed unit, current scratch owner, 마지막 stdout을 보고한다.
Epoch나 candidate grid를 조용히 줄이지 않는다.

### 18.2 After M3

Cell $c$, learner $\ell$, sieve $s$, ratio route $q$의 M3 median runtime을
$\widetilde T(c,\ell,s,q)$라 하면

$$
T_{\mathrm{single}}
=\sum_{c,\ell,s,q} n_c\widetilde T(c,\ell,s,q)
$$

로 single-worker 시간을 예측한다. 다음도 함께 보고한다.

- median, 90%, maximum runtime
- peak resident memory
- oracle Monte Carlo share
- 1, 2, 4, 8 worker measured scaling
- expected M5 wall time
- 20% operational reserve

근거: [audit plan §16](./2026-09-12-cate-experiment-audit-and-confirmatory-plan.md#16-compute-budgeting).

Measured scaling 전에는 8-core 시간이나 비용을 약속하지 않는다.

### 18.3 Count ledger

| stage | independent stochastic units | full count와 관계 |
| --- | ---: | --- |
| M0/M1 | 0 performance reps | fixture evidence |
| M2 | 6 proposed source units | smoke, 기본적으로 non-reusable |
| M3 | 420 cumulative | pilot |
| M4 | 최대 2,100 cumulative | 필요한 cell만 확대 |
| M5 | 7,800 confirmatory | synthetic, shift, controls |
| STAR | 800 | real-data source-law replications |

최종 planned confirmatory count는 M5 7,800과 STAR 800을 더한 8,600이다. 이는 model fit
수가 아니다. 현재 sequence의 M3/M4는 설계와 GO 결정에 이미 사용되므로 M5 7,800에 절대
포함하지 않는다. M5는 frozen decision record 뒤 새 disjoint seed namespace에서 시작한다.
M3/M4 결과는 pilot 또는 sensitivity evidence로 별도 보고할 수 있지만 primary M5
estimate나 confidence interval에 pooled하지 않는다.

향후 seamless internal-pilot reuse를 원하면 M3 전에 threshold, adaptation rule, stopping
rule, combination rule을 모두 동결하고 그 결합의 유효한 inference를 사전 정의해야 한다.
그 설계는 현재 sequence의 수정이 아니라 별도 protocol이다.

## 19. Execution order and authorization boundaries

실제 실행 순서는 다음과 같다.

1. Stage 0 read-only preflight
2. 별도 승인 후 current-code M0/M1 reproduction
3. 별도 승인 후 P0 implementation
4. Post-fix unit tests와 새 fingerprint M0/M1
5. 별도 승인 후 M2 six-unit smoke
6. Independent code, statistics, label review
7. 별도 승인 후 M3 420-rep pilot
8. M3 결과와 runtime을 사용해 decision record 작성
9. M4/M5 outcome 전에 threshold, multiplicity, tail cutoff freeze
10. 별도 승인 후 M4 targeted escalation
11. Terminal decision exactly one
12. `GO_CONFIRMATORY`일 때만 별도 승인 후 M5 7,800
13. Synthetic validation 완료 뒤 repaired STAR 800
14. NSW/CPS와 NSW/PSID supporting benchmark
15. WHI feasibility PASS일 때만 별도 구현·실행 계획
16. Frozen raw에서 summary, C1-C7, C1-C8 생성
17. ATE-CATE integrated table과 manuscript handoff packet
18. 별도 승인 후 manuscript, mirror, commit, push

## 20. Final go/no-go checklist

### Engineering completion

- [ ] Current fingerprint M0/M1 exists
- [ ] P0 hard tests all pass
- [ ] M2 six-unit accounting complete
- [ ] No role leakage
- [ ] No source-law or propensity identity failure
- [ ] No theorem-label contamination
- [ ] No dropped or redrawn seed
- [ ] Immutable raw and exact manifests exist

### Scientific readiness

- [ ] M3 420 units accounted
- [ ] Primary paired estimand frozen
- [ ] Interval and bootstrap convention frozen
- [ ] Multiplicity rule frozen
- [ ] Material-effect and noninferiority criteria frozen
- [ ] Tail cutoff frozen
- [ ] Runtime and storage budget approved
- [ ] M4 terminal decision is `GO_CONFIRMATORY`

### Publication readiness

- [ ] M5 7,800 units accounted
- [ ] STAR 800 units accounted or explicitly deferred
- [ ] C1-C7 and C1-C8 trace to frozen raw hashes
- [ ] ATE-CATE table does not mix unlike metrics
- [ ] Safe/forbidden claim scan passes
- [ ] Negative and failure evidence retained
- [ ] Manuscript handoff reviewed independently

이 checklist의 목적은 좋은 방향의 결과를 강제하는 것이 아니다. 올바른 estimand, 독립
평가, complete accounting을 강제한다. Joint가 RCT-only를 이기지 못한 것은 자동으로 code
failure가 아니다. Valid한 negative result도 완료된 결과다.

</details>

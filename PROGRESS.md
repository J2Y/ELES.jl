# 공개용 코드 수정 기록

## 2026-10-18 — 검토 결과 반영

대상: `ExactLayerdEkmannSolver`. 수정 전 패키지는 `../review_build/public_code_fixes_20261018/baseline/`에 보존했습니다.

| 항목 | 수정 이유 | 수정 방법 | 검증 결과 |
|---|---|---|---|
| Taylor tail 진단 | 고차 희소 계수에서 마지막 두 항이 작아도 오차가 큼 | 다항식 차수에 맞춘 계수 창으로 확대하고 오차 보증과 구분 | 고차 반례 감지, 8층 분할 해가 기준해와 약 6.4e-15 차이 |
| 이차 역추적 | 비허용 중간 후보에서 예외로 전체 중단 | 초기값 검증 유지, 중간 후보는 line search에서 거부 | 기존 중단 사례가 유한한 결과 반환, 정밀도 부족 시 성공 표시 억제 |
| 논문 재현 코드 | 수렴표·특수함수 검증·불확실성 구간 누락 | 공개 패키지 내부에 재현 경로 추가 | 독립 복사본 실행 통과, 특수함수 비교 17개 통과 |
| 노트북 실행 설정 | 등록 커널과 실행 명령 불일치 | 등록명 `eles`로 통일, 저장 출력 갱신 | 임시 커널 등록과 4개 노트북 코드·그림 재실행 통과 |

실행 명령과 측정값은 아래에 기록했습니다.

### 1차 실행 결과

- 고차 반례의 진단값이 기존 `2.66e-18`에서 `2.67e-3`으로 바뀌어 절단 문제를 감지합니다. 같은 점성함수를 8층으로 표현하면 독립 RK4 기준과 약 `6.4e-15` 차이입니다.
- 이전에 예외로 끝나던 이차 역추적은 비허용 후보를 거부하고 유한한 결과를 반환합니다. 최종 Taylor 진단이 큰 경우 `converged=false`, `refinement_required=true`를 반환합니다.
- `Pkg.test(;allow_reresolve=false)`: **16개 testsets / 149 assertions 통과**. 기존 131개와 고차 반례·탐색 후보 거부·불확실성 검증 18개를 포함합니다.
- 공개 폴더만으로 Table 2, Table 3, 95% 구간을 포함한 Figure 2를 생성했습니다. 실험의 11/10/5층과 11개 매개변수, 오차·손실이 기존 값과 일치합니다.
- mpmath 40자리 독립 검증 17개(7개 프로파일 × 1/4층 + 추가 이차근 3개)가 통과했습니다. 최대 차이는 약 `4.1e-16`입니다.
- `ChainRulesCore.ignore_derivatives`로 검증을 AD 경로 밖에 유지하고 기존 deprecation 경고를 제거했습니다.

검증 명령(패키지 루트):

```sh
julia --project=. -e 'using Pkg; Pkg.test(;allow_reresolve=false)'
julia --project=. scripts/reproduce_paper.jl
python3 reproduction/check_closed_forms.py
```

노트북의 새 코드·저장 출력 및 설치 커널 검증도 아래 최종 결과대로 완료했습니다.

### 독립 배포 폴더 확인

- 공개 폴더만 임시 디렉터리에 복사하여 Julia 재현 명령과 Python 특수함수 검증을 실행했고 둘 다 종료 코드 0으로 완료했습니다. 원고의 `p1/code`나 `figures` 폴더를 참조하지 않습니다.
- `scripts/setup_notebook.jl`은 임시 Jupyter 데이터 경로에 `eles`를 실제로 등록했고 `KernelSpecManager`에서 이름과 프로젝트 인자를 확인했습니다. README와 4개 노트북 메타데이터도 `eles`로 통일했습니다.
- `--output=DIR`로 지정한 경로에 결과가 생성되는 것을 확인했습니다. 기본 실행은 저장된 적합값을 현재 solver로 평가하며, 최적화 재실행은 `--refit` 또는 노트북 03으로 구분합니다.

### 수치·API 호환성

- `series_tail`은 마지막 두 항 대신 `max(2, degree(K)+1)`개 차수의 계수를 각 기본해별로 정규화합니다. 여전히 엄밀한 나머지항/속도 오차 상한은 아니며 차수·분할 비교가 필요합니다.
- `fit_viscosity`의 기존 `values`, `prediction`, `loss`, `optimizer` 필드를 유지하고 `diagnostics`를 추가했습니다. 최적화가 정지해도 최종 tail이 `tail_tol`을 넘으면 `converged=false`입니다. 이는 정확하지 않은 평가를 성공으로 표시하지 않기 위한 변경입니다.
- 일반 이차식과 복소근 반경 검증은 유지했습니다. 비허용 후보를 허용하도록 검사를 삭제하지 않았습니다.
- 기존 7개 프로파일 × 1/4층의 최대 RK4 차이는 `9.63e-15`, 최대 계면 잔차는 `4.44e-16`입니다. 새 정의로 계산한 이 14개 사례의 최대 tail도 `1.31e-37`로 원고의 `1.4e-37` 상한 안에 있습니다.
- Figure 2의 Monte Carlo는 Julia 고정 seed를 사용합니다. 이전 NumPy 난수열과 다르므로 구간 끝점의 작은 차이는 가능하지만, likelihood 근사·RMS noise 정의·매개변수 수는 동일합니다.

### 검증 범위

Julia 1.12.6과 기존 캐시를 사용했습니다. 새 네트워크 다운로드 설치나 다른 Julia 버전 실행은 수행하지 않았습니다. Jupyter 통신 포트는 환경 제한으로 열지 않고, 노트북 코드를 순서대로 실행해 그림을 생성하고 저장 출력을 갱신합니다. MATLAB `bvp4c`/과거 SciPy 비교는 이 수정에서 다시 실행하지 않았습니다. 외부 저장소에 게시하지 않았습니다. 최종 repository URL/DOI는 확정 후 `CITATION.cff`에 추가할 사항으로 남깁니다.

### 최종 결과

- 4개 노트북의 모든 코드 셀을 새 Julia 프로세스에서 순서대로 실행했고, 4개 그림을 생성하여 저장 출력을 갱신했습니다. 노트북 03은 실제 최적화를 다시 수행했으며 **9.6794%, 1.4479%, 1.1675%**의 점성 오차와 원고 손실값을 재현했습니다. Jacobian은 세 모델 모두 rank 11이며 Gaussian 표본은 모두 box와 layer admissibility 조건을 만족했습니다.
- 노트북 01의 과거 precompile 오류 알림을 새 실행 출력으로 교체했습니다. 4개 노트북에 실행 오류·deprecation 경고가 남아 있지 않습니다.
- 코드 변경 후 **149/149 assertions**, 독립 특수함수 비교 **17/17**, 노트북 코드·그림 **4/4**, 독립 배포 폴더 실행, `eles` 커널 등록을 확인했습니다.
- 원고 DOCX는 변경하지 않았습니다.

수정 파일은 solver·도움말, 직접 의존성/Manifest, 회귀 테스트, README/데이터 설명, 4개 노트북입니다. 새 파일은 `examples/paper_reproduction.jl`, `scripts/reproduce_paper.jl`, `reproduction/check_closed_forms.py`와 requirements, 재현용 `data/fitted_parameters.csv`, 이 진행 기록입니다.

작업공간의 상세 실행 기록은 `../review_build/public_code_fixes_20261018/`에 있습니다. 이 경로는 배포 폴더 밖의 작업 기록이며, 배포 사용자는 위 명령으로 `results/paper/`를 새로 생성할 수 있습니다.

## 2026-10-18 — 코드 변경의 원고 영향 확인

수정 전 백업 코드와 수정 후 코드를 같은 입력으로 실행하여 순방향 해, AD gradient, 역추적 최적화, Taylor 진단, 불확실성 계산을 비교합니다. 결과에 영향이 있으면 현재 ver06 원고의 해당 수치·방법·그림만 수정하고 수정본과 리뷰본을 제공합니다. 비교 결과는 이 절에 이어서 기록합니다.

### 코드 수정 전후 직접 비교 결과와 원고 반영

이번 원고 기준은 **`S1_note_ver06.docx`**입니다. 원본을 보존하고 `S1_note_ver06_corrected.docx`와 `S1_note_ver06_review.docx`로 수정 내용을 제공합니다.

| 검토 항목 | 방법 | 결과 | 원고 조치 |
|---|---|---|---|
| 순방향 해 | 수정 전 백업/수정 후 solver를 별도 Julia 프로세스에서 동일 입력으로 실행; 39개 설정·7,999개 높이 비교 | 점성, 복소 속도, 플럭스의 최대 차이 모두 **0** | 순방향 수치와 결론 유지 |
| 역추적 | 두 코드 모두 같은 관측·초기값·box·정지조건으로 새로 최적화 | 33개 최적 계수, E_K, 손실, AD/FD 결과 차이 **0**; 세 모델 모두 수렴 | Table 3의 9.68%, 1.45%, 1.17% 및 손실 유지 |
| Jacobian | 저장 적합값에서 같은 중앙차분으로 전체 행렬 비교 | 최대 차이 **0** | 민감도·국소 불확실성의 계산식 유지 |
| 진단과 탐색 | 반환 flag와 tail, 후보 거부 수를 비교 | 논문 실험에서는 비허용 후보 거부 0회, 모두 tail_tol 통과; 이차 tail은 정의 변경으로 4.25e-29 → 6.25e-29 | tail 창·line-search 거부·최종 1e-10 진단 기준 설명 추가 |
| 95% 구간 | 기존 NumPy 계산과 공개 Julia 계산 비교 | 공분산의 상대 차이는 최대 4.72e-10. 난수열/표본 생성 방식 변경으로 구간 끝점은 최대 5.55e-4 차이(구간 폭 대비 최대 약 2.09%) | 20,000개 표본·Julia MersenneTwister seed 20261017·admissibility 조건 명시; Figure 2 교체 |

모델별 재최적화 결과(수정 전과 후가 모두 동일):

| degree | E_K (%) | loss |
|---|---|---|
| 0 | 9.679449682573507 | 6.2455739057827765e-5 |
| 1 | 1.447924404225678 | 3.0314641030707297e-5 |
| 2 | 1.167504352275057 | 3.0152441569506246e-5 |

기존 적합 계수 파일과도 재최적화 결과가 모두 일치했습니다. 이번 논문 실험의 순방향/역추적 결과는 유지됩니다. 모든 입력이나 모든 최적화 문제에서 결과가 같다는 주장은 아니며, 비허용 후보가 실제로 발생하는 경우에는 새 코드의 거부 처리가 작동합니다.

원고에는 기존 검토의 남은 수식·표기 보완과 코드 영향에 따른 방법 설명만 반영했습니다. 초록, 절 구성, k0/k1/k2 표기와 핵심 결과 수치는 유지했습니다. Figure 2의 PDF/PNG는 `../figures/example3_inversion.*`에 공개 코드와 동일한 파일로 반영했습니다. Monte Carlo 끝점의 차이는 새로운 역추적 최적값의 변화가 아닙니다.

수정본/리뷰본은 원본 스타일과 수식을 유지했고, ZIP/XML 및 Pandoc 수식 읽기 검증을 통과했습니다. 리뷰의 변경 수락 텍스트가 수정본과 일치하며, 변경하지 않은 부분은 검정입니다. 실제 Microsoft Word 화면 렌더링은 수행하지 않았습니다.

상세 측정·수정 전 백업: 작업공간 `../review_build/code_manuscript_impact_20261018/`. `comparison_summary.json`, 전후 CSV, `uncertainty_comparison.json`, `word/validation.json`에 결과를 기록했습니다.

## 2026-10-18 — ver07 검토 및 리뷰본 제공

- 사용자 요청에 따라 이후 이 원고의 검토를 마치면 리뷰 DOCX도 함께 제공합니다.
- 기준은 원본 `../S1_note_ver07.docx`입니다. `../code/build_ver07_review.py`를 추가하여 서론의 선행연구 범위, 식 (9) 뒤 상수 의존성 설명, Table 2의 오타, tail 진단 창/정규화 설명, 최종 `1e-10` 기준만 반영했습니다. 수치 solver와 적합 결과는 변경하지 않았습니다.
- 생성 방법: 원본 DOCX의 스타일·수식·관계를 보존하고 실제 변경 토큰만 빨간 취소선/파랑으로 표시했습니다. 수정되지 않은 모든 글자와 수식은 검정으로 지정했습니다. 원문 문단을 복제하지 않았습니다.
- 결과: `../S1_note_ver07_review.docx` 생성, 본문 87개 블록 중 5개만 변경. 삭제 표시를 제외한 읽기는 내부 수정안과 일치하고, 추가 표시를 제외한 읽기는 원본과 정확히 일치합니다. 원본 해시도 보존했습니다.
- ZIP/XML, 스타일 참조, native 수식 260개와 Pandoc 재읽기 일치 검증을 통과했습니다. 파일 크기는 52,855바이트입니다. 실제 Microsoft Word 화면 렌더링은 수행하지 않았습니다.
- 세부 기록: `../review_build/ver07_review_20261018/changes.json`, `validation.json`. 검토 설명은 `../S1_note_ver07_check_ko.txt`에 있습니다.

### 작성자 요청 반영 — 서론 첫 문단 마지막 문장 유지

- 이유: 작성자가 해당 문장을 수정할 의사가 없다고 명시했습니다. 이후 검토와 리뷰 생성에서 수정 대상에서 제외합니다.
- 방법: `../code/build_ver07_review.py`에서 해당 수정안을 제거하고 리뷰본을 다시 생성했습니다. `../HANDOVER.md`와 검토 기록에도 유지 결정을 기록했습니다. 현재 원본에서 작성자가 이미 삭제한 Table 2 오타는 새 변경으로 표시하지 않습니다.
- 결과: 최신 원본을 기준으로 본문 변경은 3개 블록(43, 57, 65)입니다. 서론 첫 문단 전체가 원문과 같고 검정색이며 변경 표시가 없습니다. 원본 보존, 전후 읽기, ZIP/XML/스타일 및 Pandoc 수식 260개 읽기 검증을 통과했습니다. 수치 solver는 수정하지 않았습니다.

## 2026-10-18 — ver08 검토와 MATLAB 고해상도 그림 삽입

- 이유: 작성자가 ver08 검토/리뷰본과 MATLAB에서 생성한 고해상도 TIFF 그림의 문서 삽입을 요청했습니다. 기존 MATLAB 그림 스크립트는 예전 4개 모델을 사용하므로 현재 원고의 3개 모델에 맞는 생성 경로가 필요했습니다.
- 방법: `../code/export_ver08_figure_data.jl`로 현재 ELES와 저장 적합 계수를 재평가하고 공용 국소 불확실성 계산을 다시 실행했습니다. 숫자와 입력 데이터를 `../fig/data/`에 내보내고, `../code/render_ver08_figures.m`을 실제 MATLAB R2022b Update 3에서 실행했습니다. 두 그림을 1200 DPI, 8640 x 3600 RGB TIFF로 무손실 LZW 저장하고 FIG/PDF/PNG도 함께 생성했습니다.
- 출력 보완: R2024a의 로그인 대기를 확인한 뒤 기존 라이선스의 R2022b로 실행했습니다. macOS R2022b OpenGL 래스터 출력에서 글자 위치/크기 불일치를 확인하여 Painters와 명시적 여백/범례 위치로 해결했습니다. 연한 불확실성 영역은 같은 상·하한을 사용하며, 불투명한 옅은 색으로 그려 렌더러 전환을 방지했습니다.
- 코드·결과 일치: Figure 1은 기존 데이터 대비 점성 4.17e-17, 속도 4.75e-16 미만 차이이며 독립 collocation과 4.272e-14 차이입니다. Figure 2의 점성 오차 9.67945%, 1.44792%, 1.16750%와 기존 손실이 그대로 재현됐습니다. 같은 난수 설정으로 재계산한 95% 구간 끝점은 기존 공개 결과와 차이 0입니다. 최적화를 새로 실행하지 않았고 수치 solver도 변경하지 않았습니다.
- 원고 결과: ver08의 본문에는 새로 수정할 사항을 찾지 못했습니다. `../code/build_ver08_review_figures.py`로 문장·수식·표와 서론 첫 문단을 그대로 유지하고 두 그림을 기존 캡션 앞에 삽입했습니다. `../S1_note_ver08_corrected.docx`와 `../S1_note_ver08_review.docx`를 제공하며, 리뷰에서는 추가 그림만 파란 테두리로 표시했습니다. 원본은 보존했습니다.
- 검증: 두 그림을 미리보기로 확인했습니다. Word ZIP/XML/스타일/이미지 관계, 원본 본문 일치, 수식 256개 구조와 재읽기, 그림 2개 추출, TIFF 실제 크기/DPI 및 내장 바이트 일치를 검사했습니다. 실제 Microsoft Word 화면 렌더링 및 투고 페이지 수 검증은 수행하지 않았습니다.
- 기록: `../S1_note_ver08_check_ko.txt`, `../fig/README.txt`, `../fig/figure_manifest.json`, `../review_build/ver08_review_20261018/validation.json`. 공개 패키지의 solver, 최적화 로직 및 기존 수치 결과는 유지됩니다.

## 2026-10-18 — Figure 1을 매끄러운 단봉형 경계층 프로파일로 변경

- 이유: 작성자가 기존 Figure 1의 불규칙한 점성 프로파일을 더 현실적이고 매끄럽게 바꾸도록 요청했습니다.
- 방법: `data/figure1_profile.csv`와 `examples/forward_verification.jl`에 공유 입력을 정의했습니다. 비균일 10층에서 차수를 `(2,1,2,1,2,1,1,2,0,0)`으로 두고 점성과 기울기를 연결했습니다. 일반 이차 4층·일차 4층·상수 2층, 매개변수 13개 구성은 유지됩니다. 점성은 `z=0.305`에서 `0.022175`의 단일 최댓값을 갖고 상부 배경값 `0.003125`로 감소합니다. 최솟값은 `0.003`입니다.
- 계산: 기존 ELES 알고리즘으로 속도를 다시 계산하고 `../code/verify_smooth_figure1.py`의 독립 SciPy flux-form collocation과 비교했습니다. 기울기 계면 차이는 최대 `1.04e-16`, K 계면 차이는 0입니다. 최소 Taylor 반경 `1.57519`, 차수 160/240 차이 0, tail `6.72e-34`를 확인했습니다.
- 결과 변화: 입력 프로파일이 바뀌어 이전 Figure 1과 K 최대 차이 `0.0095162`, 복소 속도 최대 차이 `0.155248`입니다. solver 오류 수정이나 역산 결과 변화가 아닙니다. Figure 2의 TIFF/PNG/PDF/FIG 해시는 모두 동일합니다.
- 원고: 새 SciPy 비교 오차 `2.0617e-14`를 `2.1e-14`, Abel 상대 변동 `1.0921e-15`를 `1.1e-15`, 최소 Wronskian `10.7645`를 `10.8`로 반영했습니다. 새 프로파일 설명과 Figure 1 캡션을 보완했습니다. 나머지 본문·표와 작성자가 유지 요청한 서론 문장은 보존했습니다.
- 그림: `../code/render_ver08_figures.m`에 Figure 1만 갱신하는 실행 경로를 추가하고 실제 MATLAB R2022b에서 재실행했습니다. 1200 DPI, 8640 x 3600 LZW TIFF 및 FIG/PDF/PNG를 갱신했고 Word corrected/review에도 삽입했습니다.
- 재현성: 공개 노트북 01과 원본 생성기 `../code/build_notebooks.py`, Example 1 검증 및 과거 그림 스크립트를 공유 프로파일로 통일했습니다. 노트북 01의 코드 셀 6개와 그림을 새 Julia 프로세스에서 재실행했고, RK4 비교 `1.47e-12`, 차수 160/240 차이 0을 확인했습니다. 검증 CSV의 쉼표가 포함된 configuration을 올바르게 따옴표 처리하도록 직렬화도 보완했습니다.
- 검증: 두 Word 출력의 ZIP/XML/스타일/그림 관계, native 수식 256개와 전후 읽기 및 Pandoc 추출을 확인했습니다. 변경 본문은 Example 1과 Figure 1 캡션 두 블록뿐입니다. 수식 내 변경되지 않은 곱셈 기호 등은 검정으로 유지했고 원본 ver08 해시는 그대로입니다. 수치 solver와 역산 코드는 수정하지 않았습니다.
- 백업·증거: `../review_build/figure1_smooth_20261018/`의 `before/`, `impact.json`, `word_validation.json`, `notebook_execution.json`; `../fig/data/figure1_metrics.json`, `figure1_solver_audit.csv`. 현재 검토 기록은 `../S1_note_ver08_check_ko.txt`입니다.

## 2026-10-18 — Figure 2를 듬성듬성한 원·엑스·삼각형으로 표시

- 이유: 작성자가 Figure 2의 d=0/1/2를 각각 o/x/삼각형으로 바꾸고 점을 듬성듬성 표시하도록 요청했습니다.
- 방법: `../code/render_ver08_figures.m`의 세 모델 표시를 연결선 없는 원(o), 엑스(x), 위쪽 삼각형(^)으로 통일했습니다. 점성과 적합 속도는 기존 1001개 행에서 매 60번째 행을 선택하여 모델별 17개 점을 그립니다. 시작 행을 1/21/41로 나눠 가까운 적합값을 구별하며, 선택한 실제 좌표는 바꾸지 않습니다. 밀도는 기존 20,000개 표본/모델과 공통 60개 bin으로 같은 PDF 정규화 히스토그램 값을 계산하고, 세 번째 bin마다 하나와 최대 bin을 표시합니다. 이번 입력에서는 모델별 20개 밀도 점입니다.
- 실행: 새 `render_ver08_figures('figure2')` 경로를 실제 MATLAB R2022b Update 3에서 실행하여 Figure 1을 다시 그리지 않고 Figure 2의 TIFF/PNG/PDF/FIG만 갱신했습니다. TIFF는 RGB, 8640 x 3600, 1200 DPI, 무손실 LZW, 1,123,738바이트입니다. 참값 선, 관측값과 95% 구간은 유지했습니다.
- 수치 영향: `../fig/data/`의 모든 파일과 Figure 1의 네 출력 파일, 원본 ver08의 SHA-256이 수정 전과 같습니다. solver, 최적화, 역추적 오차 9.67945%/1.44792%/1.16750%와 불확실성 계산 자료를 변경하지 않았습니다. 그래픽 표시를 위한 행 선택이며 불확실성 표본 수를 줄이지 않았습니다.
- Word: corrected/review를 다시 생성했습니다. 이전 두 DOCX와 ZIP 내부 항목을 직접 비교하여 각각 **Figure 2 TIFF 항목 하나만 변경됨**을 확인했습니다. 본문·수식·표·캡션과 기존 변경 표시는 그대로입니다. 원본 ver08 기준의 빨간 취소선/파랑/검정과 추가 그림의 파란 테두리를 유지합니다.
- 검증: MATLAB FIG를 다시 열어 12개 모델 표시의 모양·연결선 없음·점 개수와 원본 CSV 좌표의 정확한 일치를 확인했습니다. PNG와 PDF 미리보기도 확인했습니다. Word ZIP/XML/스타일/이미지 관계, 256개 native 수식, 전후 읽기, TIFF 크기·DPI·내장 바이트 검사를 다시 수행했습니다.
- 백업·증거: `../review_build/figure2_markers_20261018/`의 `before/`, `unchanged_hashes.json`, `matlab.log`, `matlab_marker_validation.json`, `impact.json`, `word_validation.json`. 표시 방법과 재실행 명령은 `../fig/README.txt`와 `../HANDOVER.md`에도 기록했습니다.

## 2026-10-18 — Figure 2 패널 표기를 그림 안쪽으로 이동

- 이유: 작성자가 Figure 2의 (a), (b), (c)를 그림 내부에 넣도록 요청했습니다.
- 방법: `../code/render_ver08_figures.m`의 패널 표기 함수에 안쪽 배치 옵션을 추가하고 Figure 2에만 적용했습니다. 세 표기를 범례와 데이터가 없는 왼쪽 아래의 MATLAB 정규화 좌표 `(0.04, 0.12)`에 배치했습니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 Figure 2만 다시 그려 TIFF/PNG/PDF/FIG를 갱신했습니다. TIFF는 1200 DPI, 8640 x 3600, RGB/LZW, 1,122,128바이트입니다. 원·엑스·삼각형 표시와 모든 계산값 및 Figure 1은 유지했습니다. corrected/review에도 새 TIFF를 삽입했습니다.
- 확인: 미리보기를 확인했고, FIG의 세 표기 전체가 축 내부에 있으며 모든 데이터 좌표·마커·축 범위가 이전 FIG와 같은 것을 검사했습니다. 원본 ver08, Figure 1과 수치 자료의 해시를 보존했고, 각 Word ZIP에서 바뀐 항목은 Figure 2 TIFF 하나뿐입니다. 기존 Word 수식·스타일·전후 읽기 및 TIFF 내장 검증도 다시 수행했습니다.
- 기록: `../review_build/figure2_labels_inside_20261018/`의 수정 전 백업, `matlab.log`, `matlab_label_validation.json`, `impact.json`, `word_validation.json`. 배치 위치는 `../HANDOVER.md`, `../fig/README.txt`, `../S1_note_ver08_check_ko.txt`에도 반영했습니다.

## 2026-10-18 — Figure 1 패널 표기를 플롯 내부 좌상단으로 이동

- 이유: 작성자가 Figure 1의 (a), (b), (c)도 각 플롯 안쪽 좌상단에 놓도록 요청했습니다.
- 방법: `../code/render_ver08_figures.m`의 표기 함수가 위치와 수직 정렬을 받도록 바꾸고, Figure 1은 위쪽 정렬·정규화 y=0.96으로 지정했습니다. (a)는 상부 점성 곡선을 피하는 x=0.16, (b), (c)는 x=0.04입니다. Figure 2의 기존 내부 위치는 유지했습니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 Figure 1만 다시 렌더링하여 1200 DPI, 8640 x 3600 RGB/LZW TIFF(777,378바이트)와 PNG/PDF/FIG를 갱신했습니다. Word corrected/review에 같은 TIFF를 다시 넣었습니다. 수치 solver, 원고 본문·수식·표·캡션과 계산 자료는 변경하지 않았습니다.
- 확인: 미리보기에서 데이터·범례와 겹치지 않는 것을 확인하고 MATLAB에서 세 표기의 전체 영역이 축 내부임을 확인했습니다. 원본 DOCX·모든 수치 자료·Figure 2 네 파일의 해시는 수정 전과 같고, 두 Word ZIP에서 바뀐 항목은 각각 Figure 1 TIFF 하나뿐입니다. 기존 수식 256개와 Word ZIP/XML/스타일/전후 읽기/그림 내장 검증을 다시 수행했습니다.
- 기록: `../review_build/figure1_labels_inside_20261018/`의 수정 전 백업, `unchanged_hashes.json`, `matlab.log`, `impact.json`, `word_validation.json`. 현재 배치와 재실행 방법은 `../HANDOVER.md`, `../fig/README.txt`에 반영했습니다.

## 2026-10-18 — Figure 2(a)만 기존 선과 마커를 함께 표시하도록 정정

- 이유: 작성자가 마커 추가 요청은 Figure 2(a)에만 해당하며 기존 선도 함께 그리라는 뜻이었다고 명확히 했습니다. 앞선 모든 패널의 마커만 표시한 해석을 정정합니다.
- 방법: `../code/render_ver08_figures.m`에 원래 일점쇄선/파선/실선 스타일을 복구했습니다. (a)는 기존 선에 d=0 원(o), d=1 엑스(x), d=2 위쪽 삼각형(^)을 추가하고 `MarkerIndices`로 표시 간격만 줄였습니다. 상수 모델은 정확한 층 경계 점프를 유지한 계단선과 층 중앙의 원 11개, 일차·이차 모델은 전체 1001행의 선과 마커 각 17개를 사용합니다. (b)는 원래 60-bin 밀도 히스토그램 선, (c)는 원래 전체 해상도의 적합 속도 선으로 복구하고 추가 모델 마커를 제거했습니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 Figure 2만 다시 생성했습니다. 1200 DPI, 8640 x 3600 RGB/LZW TIFF는 1,123,996바이트이며 PNG/PDF/FIG도 갱신했습니다. corrected/review에 같은 TIFF를 반영했습니다. Figure 1과 두 그림의 패널 표기 위치는 유지했습니다.
- 수치 영향: 원본 ver08, Figure 1 네 파일과 모든 수치 자료의 SHA-256이 수정 전과 같습니다. solver·최적화·20,000개 불확실성 표본/모델·95% 구간과 오차 9.68%/1.45%/1.17%는 유지했습니다. 본문·수식·표·캡션은 변경하지 않았습니다.
- 확인: MATLAB FIG를 다시 열어 (a)의 선 종류·마커 종류·개수, (b)의 60개 히스토그램 bin과 원래 선, (c)의 마커 없는 1001점 모델 선을 확인했습니다. PNG 미리보기를 확인했고, Word ZIP에서 각 파일의 Figure 2 TIFF 항목 하나만 바뀐 것을 검사했습니다. 기존 수식 256개, ZIP/XML/스타일, 전후 읽기와 TIFF 내장 검증도 다시 수행했습니다.
- 기록: `../review_build/figure2_lines_markers_a_20261018/`의 수정 전 백업, `unchanged_hashes.json`, `matlab.log`, `impact.json`, `word_validation.json`. `../HANDOVER.md`, `../fig/README.txt`, `../S1_note_ver08_check_ko.txt`도 현재 지시대로 정정했습니다.

## 2026-10-18 — 두 그림의 상하좌우 바깥 여백 축소

- 이유: 작성자가 두 그림의 바깥 여백을 줄이도록 요청했습니다.
- 방법: `../code/render_ver08_figures.m`에서 실제 Painters PNG 미리보기의 내용 경계를 측정하고 상하좌우 0.04인치(약 1 mm)만 남도록 물리적 캔버스를 줄였습니다. 축과 범례를 인치 단위로 옮겨 플롯·글자 크기를 유지했습니다. TIFF/PNG/PDF를 작은 캔버스에서 새로 렌더링하고 같은 배치를 FIG에 저장합니다. 래스터 리샘플링이나 출력 뒤 이미지 잘라내기는 하지 않았습니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 두 그림을 생성했습니다. 1200 DPI RGB/LZW TIFF는 Figure 1이 8196 x 3129(746,494바이트), Figure 2가 8196 x 3211(1,118,374바이트)입니다. PNG에서 상하좌우 흰 여백은 이전 36~63픽셀에서 약 8~9픽셀로 줄었습니다. Word corrected/review도 새 종횡비로 다시 삽입했습니다.
- 유지: 수치 solver와 모든 수치 자료·원본 ver08 해시는 그대로입니다. 본문·수식·표·캡션과 Figure 1의 좌상단 표기, Figure 2의 내부 표기 및 (a)에만 선과 마커를 함께 표시하는 방식도 유지했습니다.
- 확인: 두 미리보기에서 축 제목·눈금·범례가 잘리지 않는 것을 확인했습니다. 이전 PDF와 새 PDF의 글자 내용이 같고, 모든 방향 여백이 절반 이하로 줄었음을 확인했습니다. Word 변경은 두 TIFF와 그림 높이 메타데이터에 한정됩니다. 수식 256개·전후 읽기·스타일·TIFF 내장 검증은 유지하며, 그림 크기 검사는 현재 manifest의 실제 크기와 종횡비를 확인하도록 갱신했습니다.
- 기록: `../review_build/figure_margins_20261018/`의 수정 전 백업, `before_margins.json`, `unchanged_hashes.json`, `matlab.log`, `impact.json`, `word_validation.json`. 현재 재실행 방법과 여백 기준은 `../fig/README.txt`, `../HANDOVER.md`에도 반영했습니다.

## 2026-10-18 — 그림의 변수만 이탤릭체로 통일

- 이유: 작성자가 변수명을 이탤릭체로 표시하도록 요청했습니다.
- 방법: `../code/render_ver08_figures.m`의 두 그림 축 제목·범례에서 K, z, u, v, d에만 중첩 MATLAB TeX 그룹을 적용했습니다. 설명 단어, 숫자, 괄호, 등호와 패널 표기는 정자로 두었습니다. 그룹 내부에 불필요한 공백이 생기지 않게 구성했습니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 두 그림의 TIFF/PNG/PDF/FIG를 다시 생성했습니다. 1200 DPI RGB/LZW TIFF는 Figure 1이 8201 x 3129(752,514바이트), Figure 2가 8201 x 3211(1,118,794바이트)입니다. 약 1 mm 바깥 여백과 패널 표기 위치, Figure 2(a)의 선+마커 및 나머지 패널 표시를 유지했습니다. corrected/review에도 같은 TIFF를 반영했습니다.
- 확인: PDF 글자별 폰트에서 Times-Italic은 Figure 1의 K 1개·z 3개·u/v 각 3개, Figure 2의 K 4개·z 3개·d 9개·u/v 각 1개뿐입니다. 설명 문구와 숫자는 이탤릭체에 포함되지 않습니다. 이전 PDF와 글자 내용이 같고 두 PNG의 축·범례와 가장자리 여유도 확인했습니다.
- 수치·원고 영향: solver와 모든 수치 자료·원본 ver08의 해시는 유지했습니다. Word 변경은 그림 TIFF와 종횡비에 따른 높이 정보에 한정됩니다. 원고 본문·표·수식·캡션은 유지하고, 기존 수식 256개·전후 읽기·스타일·이미지 내장 검증을 다시 수행했습니다.
- 기록: `../review_build/figure_variables_italic_20261018/`의 수정 전 백업, `matlab_compact_italics.log`, `pdf_typography_validation.json`, `impact.json`, `word_validation.json`. 표기 원칙은 `../HANDOVER.md`, `../fig/README.txt`, `../S1_note_ver08_check_ko.txt`에도 기록했습니다.

## 2026-10-18 — 최종 그림에 맞춘 ver08 본문·캡션과 리뷰본

- 이유: 작성자가 최종 그림에 맞게 원고를 수정하고 ver08_review로 보여 달라고 요청했습니다.
- 방법: `../code/build_ver08_review_figures.py`에서 원본 ver08을 유지한 채 Figure 1/2의 캡션과 Example 3 본문의 패널 참조를 보완했습니다. Figure 1(b,c)의 해석 속도/독립 collocation, (c)의 지균풍 별표를 설명했습니다. Figure 2(a)의 기존 선과 듬성듬성한 원·엑스·삼각형, (b)의 밀도 히스토그램·검정 참값 수직선, (c)의 관측점·적합선을 명시했습니다. 본문은 점성 복원을 (a), 국소 불확실성을 (a,b), 속도 비교를 (c)와 연결합니다.
- 범위: 이전 Figure 1의 C1 단봉형 프로파일 설명·검증 수치 수정은 유지했습니다. 원본 대비 변경 블록은 `[59,61,65,70,72]`이며, 이번 작업에서 이전 수정본 대비 새로 바뀐 블록은 `[61,65,70,72]`입니다. 초록·서론·표·나머지 본문과 결론은 그대로입니다. 그림 파일·그림 코드·모든 수치 자료·원본 ver08 해시는 동일하며 solver와 최적화를 실행하거나 변경하지 않았습니다.
- 결과: `../S1_note_ver08_review.docx`와 `../S1_note_ver08_corrected.docx`를 재생성했습니다. 리뷰는 실제 삭제만 빨간 취소선, 새 부분만 파랑, 공통 내용은 검정이며 문단을 복제하지 않습니다. 기존 1200 DPI TIFF 두 개는 동일한 바이트로 내장됩니다.
- 검증: 원본/수정문 두 읽기, 수식 256개와 기존 숫자 수정 두 개만 유지됨, ZIP/XML/스타일/이미지 관계, TIFF 크기·DPI·내장 바이트 및 원본·그림·자료 해시를 확인했습니다. 새 본문과 캡션이 현재 각 패널의 표시를 정확히 설명하는 것도 확인했습니다. 실제 Word 자체 화면 렌더링은 수행하지 않았습니다.
- 기록: `../review_build/ver08_figure_text_20261018/`의 수정 전 백업, `unchanged_hashes.json`, `impact.json`, `word_validation.json`; 전체 변경 사유는 `../review_build/ver08_review_20261018/figure_manuscript_changes.json`입니다. 검증 스크립트의 허용 변경 블록도 현재 범위로 갱신했습니다.

## 2026-10-18 — 패널 문자 제거와 left/middle/right 참조로 통일

- 이유: 작성자가 두 그림의 (a), (b), (c)를 없애고 캡션에서 left/middle/right로 지칭하도록 요청했습니다. 이전의 좌상단·내부 배치 지시를 대체합니다.
- 방법: `../code/render_ver08_figures.m`에서 여섯 패널 문자의 생성과 사용하지 않는 표기 함수를 제거했습니다. `../code/build_ver08_review_figures.py`의 두 캡션은 Left/Middle/Right로 구분하고 본문 Figure 2 참조도 위치로 바꿨습니다. 원본 ver08 기준 inline 리뷰 방식은 유지합니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 두 TIFF/PNG/PDF/FIG를 다시 생성했습니다. 1200 DPI RGB/LZW TIFF 크기는 Figure 1이 8201 x 3129(746,038바이트), Figure 2가 8201 x 3211(1,107,200바이트)입니다. corrected/review의 내장 그림과 캡션·본문도 갱신했습니다.
- 유지: 변수 이탤릭체, 약 1 mm 바깥 여백, Figure 2 왼쪽의 선+듬성듬성한 마커, 가운데 밀도 히스토그램과 오른쪽 적합선은 그대로입니다. 모든 수치 자료와 원본 ver08의 해시는 동일하며 solver·최적화는 변경하거나 실행하지 않았습니다.
- 검증: MATLAB FIG에 문자 표기가 없고 각 그림에 세 플롯이 남는 것을 확인했습니다. PNG를 확인하고, PDF의 글자 내용은 정확히 (a)/(b)/(c)만 줄었음을 검사했습니다. 캡션/본문의 새 읽기에 문자 패널 참조가 남지 않는 것, Word 수식 256개·ZIP/XML/스타일·리뷰 전후 읽기·1200 DPI TIFF 내장을 확인했습니다.
- 기록: `../review_build/figure_panel_names_20261018/`의 수정 전 백업, `unchanged_hashes.json`, `matlab.log`, `impact.json`, `word_validation.json`. 현재 지시는 `../HANDOVER.md`, `../fig/README.txt`, `../S1_note_ver08_check_ko.txt`에도 반영했습니다.

## 2026-10-18 — 그림 2 왼쪽 마커를 모델별 20개로 증가

- 이유: 작성자가 기존 마커가 성기게 보인다며 약 20개로 늘리도록 요청했습니다.
- 방법: `../code/render_ver08_figures.m`에서 세 모델 모두 기존 1001행 높이 자료 중 20개를 선택하도록 통일했습니다. 서로 가까운 모델은 표시 높이를 엇갈리게 합니다. 상수 모델은 선택 높이의 점을 실제 상수 구간에 삽입하여 모든 계단 점프를 보존하고, 일차·이차 모델은 전체 1001행 선에 MarkerIndices만 적용합니다. 기존 11/17/17개 마커가 20/20/20개가 되었습니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 Figure 2만 재생성했습니다. TIFF는 1200 DPI RGB/LZW, 8201 x 3189, 1,132,994바이트이며 PNG/PDF/FIG도 갱신했습니다. Word corrected/review의 내장 그림과 종횡비에 따른 높이를 갱신했습니다.
- 유지: Figure 1 네 파일과 원본 ver08, 모든 수치 자료 및 Word 생성 코드의 해시는 동일합니다. solver·최적화는 변경하지 않았습니다. 본문·수식·표·캡션과 left/middle/right 지칭, 변수 이탤릭체, 작은 바깥 여백은 유지했습니다.
- 확인: MATLAB FIG에서 모델별 MarkerIndices가 정확히 20개이며 각 마커 좌표가 기존 CSV 계산점과 1e-14 이내에서 일치함을 확인했습니다. PNG를 확인했고 PDF의 글자 내용도 동일합니다. 두 Word ZIP은 Figure 2 TIFF와 그림 높이 외 변경이 없으며 기존 수식 256개·전후 읽기·스타일·TIFF 내장 검증을 통과했습니다.
- 기록: `../review_build/figure2_20markers_20261018/`의 수정 전 백업, `unchanged_hashes.json`, `matlab.log`, `impact.json`, `word_validation.json`. 표시 수는 `../fig/figure_manifest.json`, `../HANDOVER.md`, `../fig/README.txt`에도 갱신했습니다.

## 2026-10-18 — 세 모델의 20개 마커를 같은 높이에 정렬

- 이유: 작성자가 마커를 같은 위치에 두도록 명확히 했습니다. 이전의 모델별 높이 엇갈림을 제거합니다.
- 방법: `../code/render_ver08_figures.m`에서 선택 행을 모델 반복문 밖에서 한 번 정하고, 세 모델 모두 동일한 기존 높이 `z=0.025, 0.075, ..., 0.975`를 사용하도록 바꿨습니다. 모델 간 선택 높이가 같은지 코드에서 확인합니다. 각 모델의 점성 값과 선은 그대로입니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 Figure 2만 재생성했습니다. TIFF는 1200 DPI RGB/LZW, 8201 x 3189, 1,117,386바이트이며 PNG/PDF/FIG와 corrected/review도 갱신했습니다.
- 확인: MATLAB FIG에서 세 모델이 각각 20개 마커를 가지며 마커 높이 배열이 정확히 같음을 확인했습니다. Figure 1·원본 ver08·수치 자료·Word 생성 코드의 해시는 유지했습니다. 두 Word ZIP에서 바뀐 항목은 각각 Figure 2 TIFF 하나뿐이며, 기존 본문·캡션·수식 256개와 리뷰 전후 읽기·스타일·이미지 내장 검증도 유지합니다.
- 기록: `../review_build/figure2_aligned_markers_20261018/`의 수정 전 백업, `unchanged_hashes.json`, `matlab.log`, `impact.json`, `word_validation.json`. 공통 높이는 `../fig/figure_manifest.json`, `../HANDOVER.md`, `../fig/README.txt`에도 기록했습니다.

## 2026-10-18 — 모델마다 레이어당 같은 수의 마커 표시

- 이유: 작성자가 마커의 총개수나 절대 높이 대신 레이어당 같은 개수로 표시하도록 지시했습니다.
- 방법: `../code/render_ver08_figures.m`에서 모든 모델의 각 레이어에 마커를 2개씩, 레이어 내부 1/4·3/4 위치에 배치했습니다. 상수 모델은 CSV의 실제 층 경계, 일차·이차 모델은 Julia 재현 코드의 균등 층 구성을 사용합니다. 총 마커 수는 11/10/5층에 따라 22/20/10개입니다. 상수 구간에 추가한 점은 같은 K 값을 사용하며, 일차·이차 모델은 기존 1001행 중 해당 위치를 선택하여 전체 선을 유지합니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 Figure 2만 재생성했습니다. TIFF는 1200 DPI RGB/LZW, 8201 x 3189, 1,136,560바이트이며 PNG/PDF/FIG와 corrected/review도 갱신했습니다.
- 확인: MATLAB FIG와 manifest에서 모든 레이어에 정확히 2개가 들어가고 레이어 내부 상대 위치가 1/4·3/4임을 확인했습니다. Figure 1·원본 ver08·수치 자료·Word 생성 코드의 해시는 유지했습니다. Word는 각 ZIP의 Figure 2 TIFF 한 항목만 바뀌었고 본문·캡션·수식 256개·스타일·리뷰 전후 읽기·이미지 내장 검증도 유지했습니다.
- 기록: `../review_build/figure2_per_layer_markers_20261018/`의 수정 전 백업, `unchanged_hashes.json`, `matlab.log`, `impact.json`, `word_validation.json`. 새 배치 기준은 `../HANDOVER.md`, `../fig/README.txt`, `../S1_note_ver08_check_ko.txt`에도 기록했습니다.

## 2026-10-18 — 화면상 곡선 전체의 호 길이로 마커 간격 통일

- 이유: 작성자가 '등간격'은 높이가 아니라 곡선을 따라 보이는 간격이라고 답했습니다. 전체 곡선에서의 등간격을 우선하고, 기존 총개수 22/20/10은 유지했습니다. 레이어당 2개 조건은 전곡선 등간격과 일반적으로 양립하지 않으므로 레이어별 개수는 달라진다고 알렸습니다.
- 방법: `../code/render_ver08_figures.m`에 물리적 축 가로·세로 크기와 축 범위를 반영하는 호 길이 배치를 추가했습니다. 전체 표시 선의 누적 길이를 구해 등분하고 양 끝에는 반 간격을 남겼습니다. 마커 위치는 기존 선분 위에 삽입하므로 원래 선의 꺾임·모양과 일차/이차의 원래 1001개 꼭짓점이 유지됩니다. 상수 모델은 수평 계단 연결선도 포함합니다.
- 표기의 의미: 마커는 모델 선을 구분하는 그림 기호이며 추가 계산/관측 표본이 아닙니다. 계단 연결선의 마커가 중간 점성값을 뜻하지 않도록 `../code/build_ver08_review_figures.py`의 Figure 2 캡션에 모델 구분·연결선 포함·화면 좌표 기준을 명시했습니다. 그 외 본문·표·수치·캡션은 유지했습니다.
- 실행·결과: 실제 MATLAB R2022b Update 3에서 Figure 2만 재생성했습니다. 1200 DPI RGB/LZW TIFF는 8201 x 3189, 1,135,034바이트이며 PNG/PDF/FIG와 Word corrected/review도 갱신했습니다. 실제 축에서 마커 간 호 길이는 모델별 약 0.205170/0.177219/0.354384인치입니다.
- 검증: 완성 FIG를 다시 열어 실제 축 크기와 선분을 기준으로 측정했습니다. 등분 목표 위치와의 최대 차이는 2.22e-16인치이며 22/20/10개가 맞습니다. PNG를 확인했고 PDF의 글자 내용은 동일합니다. Figure 1·원본 ver08·모든 수치 자료 해시는 유지했습니다. Word 변경은 Figure 2 TIFF와 캡션의 마커 설명 및 해당 대체텍스트뿐입니다. native 수식 256개·스타일·리뷰 전후 읽기·TIFF 내장 검증을 통과했습니다.
- 기록: `../review_build/figure2_visual_spacing_20261018/`의 수정 전 백업, `unchanged_hashes.json`, `matlab.log`, `matlab_arc_validation.json`, `impact.json`, `word_validation.json`. 현재 원칙은 `../HANDOVER.md`, `../fig/README.txt`, `../S1_note_ver08_check_ko.txt`에도 반영했습니다.

## 2026-10-19 — README에 논문 제목과 AML 투고 상태 표시

- 이유: 작성자가 README에 논문 제목과 Applied Mathematics Letters에 투고한 원고라는 설명을 넣도록 요청했습니다.
- 방법: 최신 `../S1_note_ver09.docx`의 제목 `A differentiable exact solver for layered Ekman flows with piecewise-continuous eddy viscosity`를 확인하여 공개용 `README.md`와 원고 폴더의 `../README.md`에 반영했습니다. 작성자의 지시에 따라 상태를 `Submitted to Applied Mathematics Letters (AML)`로 표시했습니다.
- 결과: 두 README에서 같은 제목과 투고 저널을 확인할 수 있습니다. solver, 수치 자료, 그림과 원고 파일은 수정하지 않았으며 기존 계산 결과에 영향이 없습니다.
- 확인: 저장된 두 README의 제목과 투고 상태를 다시 읽어 확인했습니다. 문서 설명만 바꿨으므로 수치 실험은 다시 실행하지 않았습니다.

## 2026-10-18 — 현재 ver08의 AML 출판 분량 확인

- 이유: 작성자가 현재 논문이 Applied Mathematics Letters의 출판 형식에서 6페이지에 들어가는지 확인하도록 요청했습니다. 과거 HANDOVER의 2단 5p 가정과 오래된 원고의 페이지 수만으로 판단할 수 없어 현재 Word 파일에서 다시 조판했습니다.
- 방법: `../review_build/ver08_aml_length_20261018/build_length_check.py`를 추가했습니다. 현재 corrected Word의 모든 80개 Pandoc 블록과 수식 256개를 보존하고, 수식용 Word 표를 실제 equation으로, 표 3개와 그림 2개를 float로 변환했습니다. 본문과 캡션은 줄이지 않았습니다. 최근 AML 게재본의 실제 PDF에서 페이지 크기·본문 폭·서체 크기·줄간격을 측정했고, 공개 OFL 글꼴 Charis SIL을 감사 폴더에만 내려받아 사용했습니다. 첫 페이지의 저널 머리말과 출판 메타데이터 공간도 반영했습니다.
- 결과: 일반 elsarticle 3p/10pt 템플릿은 8페이지입니다. 실제 게재본에서 측정한 본문/참고문헌 크기로 조판하고 출판 부가 공간을 남기면, 참고문헌 학술지명을 풀어 쓴 버전은 7페이지이고 학술지명만 표준 약어로 바꾸면 **6페이지**입니다. 저자·논문 제목·연도·권호·쪽수·표시된 DOI는 모두 유지했습니다. 전달 파일은 `../S1_note_ver08_AML_layout_check.pdf`이며 출판사의 공식 proof는 아닙니다. 마지막 쪽 여유가 작습니다.
- 확인: 모든 수식 payload 256개와 번호식 13개, 표 3개, 그림 2개, 참고문헌 8개를 확인했습니다. 학술지명 약어를 되돌리면 생성한 본문 텍스트가 같고, 전달 PDF에는 overfull box와 누락 글자가 없습니다. 페이지 2/6의 표·수식·그림·참고문헌 배치를 확인했습니다. 원본/수정본/리뷰본 Word, 기존 Manuscript_AML.tex와 그림 PDF의 해시는 그대로입니다.
- 수치·원고 영향: solver와 재현 실험 코드는 수정하지 않았고 계산을 다시 실행하지 않았습니다. 기존 수치 결과와 Word 본문·수식·그림은 변하지 않습니다. 이번 변경은 분량 확인용 조판 코드·PDF 및 감사 기록에 한정됩니다.
- 기록: `../review_build/ver08_aml_length_20261018/`의 `published_geometry.json`, `published_reference_sources.json`, `input_hashes.json`, `conversion_manifest.json`, `validation.json`, 각 TEX/PDF/로그와 미리보기. `../HANDOVER.md` 맨 위에 현재 결과와 과거 5p 가정의 한계를 기록했습니다.

## 2026-10-19 — ver09 정규화 좌표의 Taylor 계수 hat 확인

- 이유: 작성자가 식 (11)의 정규화 좌표 계수를 표 1의 원좌표 계수와 구분하도록 hat을 요청했습니다.
- 방법: `../review_build/ver09_hat_coefficients_20261019/build_hat_review.py`를 추가하여 식 (11) 위의 다항식, 식 (11), 바로 뒤의 계수 목록·중심값 정의·일반화된 전개와 재귀식에 해당하는 13개 기호를 검사합니다. hat이 없는 계수만 native Word 수식 accent로 보완하고, 리뷰는 실제 추가 변경에 한해 기존 k를 빨간 취소선, 새 hat k를 파랑으로 표시합니다. 첨자와 공통 부분은 유지하며 표 1과 층 번호에는 hat을 추가하지 않습니다.
- 결과: 작업 중 작성자 원본이 외부에서 갱신되었습니다. 최종 기준 `../S1_note_ver09.docx`(SHA256 `7f5200a9720ef2a0b87e369c3eb7805bf502b3da654408e4b18fc3bb5f4e3178`)에는 13곳 모두 이미 hat이 있어 추가 수정은 0곳입니다. 최신 원본을 보존하고 같은 본문으로 `../S1_note_ver09_corrected.docx`, `../S1_note_ver09_review.docx`를 제공했습니다. 리뷰에 기존 반영 사항을 새 변경처럼 표시하지 않았습니다.
- 검증: 두 파일의 native 수식 255개와 Pandoc으로 읽은 hat k 13개, 그림 2개, ZIP/XML·스타일 참조를 확인했습니다. 원본의 모든 본문 블록과 수식·그림·스타일 등 각 ZIP 항목의 내용이 그대로이고, 리뷰의 전후 읽기가 원본·corrected와 일치합니다. 원본 해시를 다시 확인했습니다. 실제 Microsoft Word 화면 렌더링은 수행하지 않았습니다.
- 수치 영향: 수치 solver, 재현 실험 코드와 계산 자료는 변경하지 않았습니다. 표기 확인·Word 패키징용 코드만 추가하여 수치 실험을 다시 실행하지 않았으며 논문의 결과 수치는 그대로입니다.
- 기록: `../review_build/ver09_hat_coefficients_20261019/validation.json`, 입력 스냅샷과 Pandoc 재읽기 TEX에 있습니다.

## 2026-10-19 — J2Y/ELES.jl GitHub 초기 배포 준비 및 검증

- 이유: 작성자가 공개용 패키지를 검토한 뒤 `https://github.com/J2Y/ELES.jl`에 푸시하도록 명시적으로 요청했습니다.
- 수정 방법: README에 실제 저장소 복제·설치 명령과 검증 버전 Julia 1.12.6을 추가했습니다. CITATION의 이전 폴더명 표기를 정리하고 `repository-code`를 실제 저장소 URL로 설정했습니다. 노트북 01의 저장된 실행 알림에서 로컬 OneDrive 경로 한 곳만 다른 노트북과 같은 `<package>`로 치환했습니다. 노트북의 모든 코드·수치·그림·메타데이터는 그대로입니다.
- 배포 범위: 소스, 테스트, 재현 스크립트, CSV 입력 자료, 실행 출력이 저장된 노트북 4개, README·인용 정보·MIT 라이선스·이 진행 기록을 포함한 25개 파일입니다. 재생성되는 `results/`와 로컬 검토 산출물은 포함하지 않습니다. 비밀값·로컬 절대경로 패턴, 외부 파일 의존, symlink 및 대용량 파일을 점검했고 배포를 막을 문제는 발견하지 못했습니다.
- 실행·결과: 별도 복사본에서 Julia 1.12.6으로 **16개 testsets / 149 assertions 통과**, `scripts/reproduce_paper.jl` 기본 실행 성공, 독립 40자리 해석해 검증 **17/17 통과**(최대 차이 `4.07e-16`)를 확인했습니다. 새 표 2·3, 카탈로그 지표·속도, 국소 구간 CSV 5개는 기존 결과와 바이트 단위로 같습니다. 점성 오차는 **9.6794497%, 1.4479244%, 1.1675044%**로 유지됩니다.
- 수치·원고 영향: 수치 solver와 재현 실험 코드는 수정하지 않았고 원본 결과 8개 파일의 해시가 유지됩니다. 따라서 이번 배포 정리로 원고 수치나 그림을 수정할 필요는 없습니다. 이번에는 전체 노트북 재실행이나 `--refit` 재최적화를 수행하지 않았습니다.
- 실행 환경: 최초 전체 자동 precompile은 선택적 ZMQ의 소켓 생성이 sandbox에서 제한되어 실패했습니다. `JULIA_PKG_PRECOMPILE_AUTO=0`으로 불필요한 전체 사전 컴파일을 끈 뒤 수치 테스트와 재현을 완료했습니다. 의존성은 기존 캐시를 사용했습니다.
- Git 방법·근거: 기존 SSH 인증으로 빈 원격 저장소를 확인하고 별도 임시 clone의 `main`에서 첫 커밋을 준비했습니다. 기존 작업 폴더를 덮어쓰지 않으며 강제 푸시는 사용하지 않습니다. 게시 커밋과 원격 확인 결과는 작업공간 `../review_build/github_release_20261019/release_state.json`, 검증 로그는 같은 폴더의 `validation/`에 기록합니다.

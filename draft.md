# Nhánh 1 — Tiền xử lý dữ liệu biểu hiện gene (RNA-seq)

Tài liệu mô tả **toàn bộ quá trình tiền xử lý của nhánh dữ liệu gene expression (RNA-seq)**: định nghĩa format cuối, từng bước xử lý kèm và số mẫu còn lại của từng cohort sau mỗi bước.

---

## 0. Định nghĩa nhánh và phạm vi tài liệu

| | |
|---|---|
| **Tên nhánh** | `BRANCH_EXPR`|
| **Tiêu chí** | `qc_has_rnaseq = True` — có profile biểu hiện toàn hệ gene tại **thời điểm chẩn đoán**, một mẫu đại diện cho một bệnh nhân |
| **Phạm vi tài liệu này** | phần **RNA-seq sau khi merge**: 1.252 bệnh nhân / 5 cohort |
| **Chưa nằm trong tài liệu này** | 721 bệnh nhân microarray (GSE10358 304, GSE12417 242, GSE71014 104, GSE5122 58, GSE106291 13)|

---

## 1. FORMAT

### 1.1. Thành phần A — ma trận biểu hiện

| Mục | Đặc tả |
|---|---|
| File | `AML_CPM_pyComBat.parquet` |
| Nội dung | ma trận biểu hiện gene **đã chuẩn hoá và đã hiệu chỉnh batch**, dùng trực tiếp làm feature matrix |
| Biểu diễn ma trận | **gene × bệnh nhân** (hàng = gene, cột = bệnh nhân) |
| Khoá hàng | `gene_id` = Ensembl gene ID (vd `ENSG00000000419`) |
| Khoá cột | `merged_id` dạng `<Cohort>\|<sample_id>` (vd `BeatAML\|BA2409`) |
| Kích thước | **12.058 gene × 1.252 bệnh nhân** |
| Đơn vị | log₂(1 + CPM) |
| Kiểu dữ liệu | `float32`; **0 ô NaN** |
| Số thuộc tính (feature) | **12.058**, trong đó `protein_coding` **11.671** (96,8 %), còn lại 387 gene thuộc các biotype khác (pseudogene 335, IG gene 47, TR gene 4, lncRNA 1) |
| Bảng tra gene | `AML_gene_annotation.csv` — `gene_id`, `symbol`, `biotype` |

### 1.2. Thành phần B — bảng label, **14 thuộc tính**

File `AML_labels_aligned_CPM_pyComBat.csv`, 1.252 dòng × 14 cột. Một dòng = một bệnh nhân.

| # | Cột | Nội dung | Kiểu / miền giá trị | Độ phủ |
|---|---|---|---|---|
| 1 | `merged_id` | khoá nối với cột ma trận | string `<Cohort>\|<sample_id>` | 1.252 (bắt buộc) |
| 2 | `patient_id` | ID bệnh nhân toàn cục | string `<COHORT>:<id gốc>` | 1.252 |
| 3 | `cohort` | tên cohort | category, 5 giá trị | 1.252 |
| 4 | `sample_id` | ID mẫu RNA-seq trong file nguồn | string | 1.252 |
| 5 | `timepoint` | thời điểm lấy mẫu | category — toàn bộ = `Diagnosis` | 1.252 |
| 6 | `os_time_days` | thời gian sống còn toàn bộ | numeric, **ngày** (trung vị 350) | 1.243 |
| 7 | `os_event` | biến cố tử vong | 1 = tử vong (760) / 0 = censored / NA | 1.243 |
| 8 | `cr` | đạt CR sau induction | 1 = CR/CRi / 0 = không đạt / NA | 953 |
| 9 | `eln2017` | phân nhóm nguy cơ ELN 2017 | Favorable 310 / Intermediate 259 / Adverse 409 / Unknown 274 | 978 (ngoài `Unknown`) |
| 10 | `age_years` | tuổi lúc chẩn đoán | numeric, trung vị 61 | 1.250 |
| 11 | `sex` | giới tính | F 606 / M 646 | 1.252 |
| 12 | `source_cohort` | cohort nguồn — **biến batch cấp 1** | category, 5 giá trị | 1.252 |
| 13 | `wave` | đợt giải trình tự (chỉ BeatAML) | `Waves1+2` / `Waves3+4` / `Both` / NA | 456 |
| 14 | `batch_wave` | biến batch dùng cho pyComBat | `BeatAML_w12` 301 / `ClinSeq` 315 / `GSE146173` 246 / `BeatAML_w34` 155 / `TCGA-LAML` 151 / `FIMM` 84 | 1.252 |

**Quy ước giá trị thiếu:** `NA` = không rõ / cohort không có trường đó.

---

## 2. CÁC BƯỚC XỬ LÝ

### B0 · Thu thập nguồn

| Cohort | n bệnh nhân cuối | File nguồn chính |
|---|---|---|
| BeatAML 2.0 | 456 | `BeatAML2_full.RData` (`rnaCounts`, raw counts) |
| ClinSeq | 315 | `SSFAML_Clinseq_Salmon_GeneExpression_geneCount.csv` |
| GSE146173 | 246 | `GSE146173_GEO_raw_gene_counts.txt.gz` |
| TCGA-LAML | 151 | `laml.rnaseq.179_v1.0_gaf2.0_read_count_matrix.txt` |
| FIMM-AML | 84 | `File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv` |

### B1 · Trích xuất từng cohort về một schema chung + lọc mẫu

Quy tắc: chỉ giữ mẫu **thời điểm chẩn đoán**, **một mẫu / một bệnh nhân** (giữ mẫu Bone Marrow Aspirate nếu trùng) đưa về `<COHORT>:<id gốc>`, 10 trường schema (`patient_id`, `cohort`, `sample_id`, `timepoint`, `os_time_days`, `os_event`, `cr`, `eln2017`, `age_years`, `sex`).

| Cohort | mẫu trong file nguồn | bộ lọc đã áp | còn lại |
|---|---|---|---|
| BeatAML | **698** mẫu RNA-seq (632 BN) | `diseaseStageAtSpecimenCollection == "Initial Diagnosis"` ∧ `isRelapse == FALSE` ∧ có mặt trong `rnaCounts` → 458 mẫu / 456 BN; 2 BN có 2 mẫu chẩn đoán → giữ mẫu Bone Marrow Aspirate| **456** |
| ClinSeq | **374** cột mẫu | khớp `patient_id` trong `ClinSeq_Clinics` ∧ loại BN có cờ `excluded` ∧ khử trùng cột | **315** |
| GSE146173 | **246** mẫu | lọc `os_time_days` không rỗng (không loại ai) | **246** |
| TCGA-LAML | **151** mẫu| chỉ `Primary Blood Derived Cancer` (máu ngoại vi + tuỷ) | **151** |
| FIMM | **167** mẫu (163 AML + 4 người hiến khoẻ mạnh) | bỏ 4 mẫu `Healthy`; `Disease status == "Diagnosis"`; 1 mẫu/BN | **84** |
| **Tổng** | **1.636** | | **1.252** |


### B2 · Hợp nhất định danh gene và lấy giao

- Cắt hậu tố phiên bản ENSG (`ENSG00000001.5 → ENSG00000001`) 
- Lấy **giao** của 5 tập gene, không hợp — gene chỉ có ở một phần cohort là nguồn phá hoại chính của bước hiệu chỉnh batch (xem B4).

| Cohort | gene đưa vào | gene còn lại nếu bỏ cohort này | gene riêng cohort này làm mất |
|---|---|---|---|
| ClinSeq | 35.547 | 49.041 | **16.906** |
| BeatAML | 56.828 | 33.129 | 994 |
| GSE146173 | 57.033 | 33.574 | 1.439 |
| TCGA-LAML | 60.616 | 32.135 | 0 |
| FIMM | 60.619 | 32.135 | 0 |
| **Giao 5 cohort** | | | **32.135** (18.730 `protein_coding`) |

→ Output: `AML_merged_counts.parquet` (32.135 × 1.252, **raw counts, chưa chuẩn hoá**)

### B3 · Chuẩn hoá thuật ngữ (mapping)

| Nhãn | Quy tắc  |
|---|---|
| `cr` BeatAML | `responseToInductionTx`: `Complete Response`, `Complete Response i` → 1; `Refractory` → 0; còn lại NA |
| `cr` GSE146173 | `CR`/`CRi` → 1; `refractory_disease`, `death_with_residual_disease` → 0; `death_in_aplasia`, `death_from_indeterminate cause` → **NA** (chết trước khi đánh giá được đáp ứng ≠ không đạt CR) |
| `cr` FIMM | `First_cr_pt`có`1` (116 BN) và NA (70 BN) |
| `cr` TCGA-LAML | **không tồn tại** `cr = NA` |
| `os_event` ClinSeq | `status` là 0/1 không nhãn; **1 = tử vong** |
| `os_time_days` TCGA | `days_to_death` nếu tử vong, nếu không `days_to_last_followup`, dự phòng `days_to_last_known_alive` |
| `eln2017` | `FavorableOrIntermediate` / `IntermediateOrAdverse` → **Unknown** (không gán vào một trong hai nhóm) |
| đơn vị OS | quy về **ngày** ở mọi cohort |

→ Output: `AML_merged_labels.csv` (1.252 × 12)

### B4 · Lọc gene — **trong từng cohort**, trước khi hiệu chỉnh batch

**Quy tắc:** giữ gene có `CPM > 1` ở **≥ 25 % số mẫu của *mỗi* cohort**. Kết quả: **32.135 → 12.058 gene**.

| Nhánh | Bộ lọc | Số gene | Kết quả |
|---|---|---|---|
| A | không lọc | 32.135 | **ma trận hỏng**: 76,0 % ô là NaN, 5/6 batch mất sạch, 951/1.252 BN không còn ô dùng được — pyComBat **không báo lỗi** |
| **B** | **lọc trong từng cohort** | **12.058** | kNN 0,032; đối chứng dương 0,990; LOCO CR AUC 0,700 |

Cơ chế: gene có phương sai bằng 0 trong ít nhất một cohort ⇒ δ(g,b) không xác định; empirical Bayes gộp thông tin mọi gene trong cùng batch để ước lượng prior nên chỉ một nhóm gene không hữu hạn là hỏng prior của cả batch.

### B5 · Chuẩn hoá — CPM trong từng cohort, rồi log₂(1 + CPM)

### B6 · Hiệu chỉnh batch — pyComBat với reference batch

| Tham số | Giá trị |
|---|---|
| Công cụ | pyComBat (`inmoose` 0.9.1) — ComBat empirical-Bayes, hiệu ứng batch cộng γ(g,b) và nhân δ(g,b) |
| `batch` | `batch_wave`: **BeatAML tách thành 2 batch theo đợt giải trình tự** — `BeatAML_w12` (n = 301, gồm 12 BN "Both") và `BeatAML_w34` (n = 155); 4 cohort còn lại mỗi cohort một batch |
| `ref_batch` | `"BeatAML_w12"`|
| `covar_mod` | `None` — **không đưa outcome vào mô hình hiệu chỉnh** |
| `par_prior` | `True`; dữ liệu vào ở thang log₂ |

### B7 · Sắp nhãn theo đúng thứ tự cột

---

## 3. TỪNG COHORT GỐC — xử lý, mapping, làm sạch và số mẫu giữ lại

### 3.1. Bảng tổng hợp in/out

| Cohort | mẫu nguồn | sau lọc timepoint/QC | **vào ma trận** | gene đưa vào | có OS | có CR | có ELN2017 | tử vong |
|---|---|---|---|---|---|---|---|---|
| BeatAML 2.0 | 698 | 458 (456 BN) | **456** | 56.828 | 456 | 332 | 441 | 249 |
| ClinSeq | 374 | 315 | **315** | 35.547 | 315 | 308 | 299 | 225 |
| GSE146173 | 246 | 246 | **246** | 57.033 | 246 | 229 | 238 | 149 |
| TCGA-LAML | 151 | 151 | **151** | 60.616 | 142 | **0** | **0** | 97 |
| FIMM-AML | 167 | 84 | **84** | 60.619 | 84 | 84 |**0** | 40 |
| **Tổng** | **1.636** | **1.254** | **1.252** | giao 32.135 → lọc 12.058 | **1.243** | **953**| **978** | **760** |

### 3.2. BeatAML 2.0 — 698 → 456

- Nguồn `BeatAML2_full.RData`
- Lọc: `Initial Diagnosis` ∧ `isRelapse = FALSE` ∧ có trong `rnaCounts` → 458 mẫu / 456 BN; 2 BN có 2 mẫu chẩn đoán → giữ mẫu Bone Marrow Aspirate.
- Batch: tách `BeatAML_w12` 301 / `BeatAML_w34` 155 theo cột `cohort` của `beataml_wv1to4_clinical.xlsx`, map qua `dbgap_rnaseq_sample`.

### 3.3. ClinSeq — 374 → 315

- Counts là **Salmon estimated counts (số thực)** ⇒ `round()`.
- Loại BN có cờ `excluded` trong `ClinSeq_Status`; khử cột trùng; chỉ giữ cột khớp `patient_id` của bảng clinical. 
- `os_event` kiểm tra status 0/1 bằng trung vị `follow_up_days` (198 vs 1.344 ngày).

### 3.4. GSE146173 — 246 → 246

- Mapping `response_to_induction_treatment` 6 trạng thái → `cr` 3 trạng thái  B3); 17 BN rơi vào NA.
- ENSG có hậu tố phiên bản → cắt trước khi giao gene.

### 3.5. TCGA-LAML — 151 → 151

- Dùng **GDC STAR-Counts (GENCODE v36)**, cột `unstranded`; bỏ các dòng `N_unmapped`, `N `
- Chỉ lấy mẫu `Primary Blood Derived Cancer` (máu ngoại vi + tuỷ).
- Không có nhãn CR và không có ELN2017 ⇒ chỉ dùng được cho endpoint OS (142/151 BN có OS).

### 3.6. FIMM-AML — 167 → 84

- Ma trận gốc chứa **4 người khoẻ mạnh** —  bỏ.
- **Ba lớp ID**: clinical khoá `Patient_ID` (`AML_019`), annotation có cả `Patient_ID` lẫn `Sample_ID`, ma trận khoá `Sample_ID` (`AML_084_01`). Join trực tiếp `Patient_ID` khớp **0 dòng** ⇒ phải join qua annotation và so bằng phần chữ số của ID.
- Chỉ giữ `Disease status == "Diagnosis"`, 1 mẫu/BN → 84.

---

## 4.  ĐÁNH GIÁ MA TRẬN CUỐI

| Nhóm | Chỉ số | Hướng tốt | Trước hiệu chỉnh | **Sau (CPM + pyComBat, ref `BeatAML_w12`)** |
|---|---|---|---|---|
| **Batch** | PC1 variance do cohort | thấp | 0,8778 | **0,0000** |
| | PC1–5 variance do cohort | thấp | 0,5850 | **0,0000** |
| | kNN batch excess (k = 30) | thấp | 0,5975 | **0,0319** |
| | Silhouette theo cohort | thấp | +0,1586 | **−0,0245** |
| | Chênh median *ACTB* giữa cohort | thấp | 1,8589 log₂ | **0,0259 log₂** |
| | Chênh median toàn ma trận | thấp | 1,5477 log₂ | **0,0251 log₂** |
| **Sinh học** | Giữ cấu trúc nội bộ cohort (TB / tệ nhất) | cao | 1,0000 | **0,9932 / 0,9891** |
|  **Đối chứng dương**| — AUC giới tính (LOCO) | cao | 0,9898 | **0,9897** |
| **Đối chứng âm** | Silhouette khi xáo nhãn cohort | ≈ 0 | −0,0183 | **−0,0210** |
| **Downstream** | LOCO CR AUC (4 cohort, 934 BN) | cao | 0,6952 ± 0,0292 | **0,6999 ± 0,0327** |
| | LOCO OS C-index (5 cohort, 1.243 BN) | cao | 0,6437 ± 0,0195 | **0,6621 ± 0,0216** |
| **Ổn định** | % ô bị đổi khi thêm cohort thứ 5 | thấp | — | **0,000 %** |

# Nhánh 2 — Tiền xử lý dữ liệu Clinical + Mutation (BRANCH_MUT)

Mô tả quá trình tiền xử lý của nhánh dữ liệu Clinical + Mutation: định nghĩa format cuối, từng bước xử lý và số bệnh nhân còn lại của từng cohort sau mỗi bước.

---

## 0. Định nghĩa nhánh và phạm vi

| | |
|---|---|
| **Tên nhánh** | `BRANCH_MUT`|
| **Tiêu chí vào nhánh** | `qc_has_mutation = True`|
| **Ràng buộc mẫu** | mẫu tại **thời điểm chẩn đoán**, 1 dòng = 1 bệnh nhân, có **ít nhất một outcome** (OS hoặc CR) |
| **Cỡ mẫu thực tế** | **5.078 bệnh nhân / 7 cohort**|

---

## 1. FORMAT

Format cuối là **bộ ba file**, nối bằng khoá `merged_id = <source_cohort>__<patient_id>`.

### 1.1. Thành phần A — `aml_core.csv` (5.078 × 58)

| Mục | Đặc tả |
|---|---|
| Nội dung | bảng **1 dòng = 1 bệnh nhân**, định danh + clinical + outcome + cờ QC |
| Khoá | `merged_id`; bắt buộc mọi dòng: `patient_id`, `source_cohort` |
| Số thuộc tính | **58 trường core** — `identifier` 8 · `clin_` 26 · `out_` 17 · `qc_` 6 · `reference` 1 |
| `timepoint` | **`Diagnosis` cho cả 5.078 dòng** |

Độ phủ thực tế của các trường chính (số dòng không rỗng / 5.078):

| Feature | Nội dung | Độ phủ |
|---|---|---|
| `clin_age_years` | tuổi lúc chẩn đoán | 5.074 |
| `clin_sex` | giới tính | 5.078 |
| `clin_wbc` · `clin_hb` · `clin_platelets` | bạch cầu · hemoglobin · tiểu cầu | 4.870 · 4.763 · 4.807 |
| `clin_bm_blast_pct` | % blast tuỷ | 4.530 |
| `clin_ecog` | toàn trạng ECOG | 4.029 |
| `clin_fab` | phân loại FAB | 819 (chỉ BeatAML, TCGA-LAML, GSE146173, FIMM) |
| `clin_karyotype_raw` | chuỗi karyotype ISCN gốc | 1.351 (UK-NCRI và AMLSG **không phát hành** chuỗi ISCN) |
| `clin_eln2017` | ELN 2017 | 4.820 (Favorable 1.833 / Intermediate 1.178 / Adverse 1.809) |
| `out_os_days` + `out_os_event` | sống còn toàn bộ (ngày) | 5.064 · 5.066 (tử vong 3.173) |
| `out_cr_binary_strict` | **nhãn CR chính** | **4.592** (CR 3.222 / không CR 1.370) |
| `out_cr_binary` | CR định nghĩa rộng | 2.433 — thấp vì UK-NCRI không tách được `early_death` (CR-02) |
| `out_cr_response` | đáp ứng nhiều mức | CR 3.165 · CRi 57 · RD 671 · `RD_or_early_death` 557 · early_death 142 · no_treatment 53 · PR 7 · not_evaluable 2 · unknown 131 · not_reported 293 |
| `qc_split` | train/test 80/20 phân tầng theo cohort, seed 42 | train 4.062 / test 1.016 |

### 1.2. Thành phần B — `aml_genomic_wide.parquet` (5.078 × 421)

| Mục | Đặc tả |
|---|---|
| Nội dung | ma trận **bệnh nhân × đặc trưng di truyền**, mỗi trường giá trị đi kèm một cột mask |
| Cấu trúc cột | 1 cột khoá `merged_id` + **210 trường giá trị** + **210 cột `assayed_<field>`** |
| Giá trị | `1` = Positive · `0` = Negative (đã xét nghiệm, âm tính) · `NaN` = **không nằm trong panel của cohort** |
| Mask | `assayed_X = 1` ⟺ ô `X` có giá trị quan sát được. Toàn ma trận **57,4 % ô được xét nghiệm** |

**210 trường theo khối:**

| Tiền tố | Số trường | Nội dung | Nguồn danh sách |
|---|---|---|---|
| `mut_<GENE>` | **119** | gene có đột biến somatic hay không (mức gene) | **universe cố định** = 119 gene `mut_` của codebook |
| `var_<GENE>_<VARIANT>` | **14** | biến thể cụ thể mà ELN 2017 phân biệt | FLT3 ITD/TKD/other, CEBPA bi/mono, IDH1 R132, IDH2 R140/R172, NRAS G12_13/Q61_62, U2AF1 S34/Q157, KMT2A PTD |
| `fus_<A>_<B>` | **14** | chuyển đoạn / gene fusion | RUNX1-RUNX1T1, CBFB-MYH11, PML-RARA, KMT2A-r (+MLLT3/MLLT1/ELL), DEK-NUP214, GATA2-MECOM, BCR-ABL1, NUP98-r (+NSD1), KAT6A-CREBBP, RBM15-MRTFA |
| `cyto_*` | **63** | bất thường nhiễm sắc thể | 62 cột của UK-NCRI/AMLSG + `cyto_normal` |

 BeatAML và TCGA-LAML là WES/WGS nên gọi được đột biến ở 2.477 và 1.449 gene. Nếu lấy hợp tất cả thì ma trận phồng lên 3.638 cột mà 3.479 cột chỉ có ở đúng một cohort.

### 1.3. Thành phần C — `aml_branch_membership.csv` (5.078 × 9)

Cột: `merged_id`, `source_cohort`, `qc_has_rnaseq`, `qc_has_mutation`, `qc_has_cytogenetics`,
`has_outcome`, `branch_MUT`, `branch_EXPR`, `branch_BOTH`, `qc_split`.
Kết quả: **branch_MUT 5.078** · branch_EXPR 1.068 · branch_BOTH 1.068.
`branch_EXPR` ở đây chỉ đếm bệnh nhân **trong 7 cohort này** có RNA-seq; các cohort chỉ có expression
(GSE12417, GSE71014, GSE5122, GSE37642…) chưa gộp vào.

### 1.4. Bốn tầng feature theo độ phủ panel

| Tầng | Số trường | `mut_` | `var_` | `fus_` | `cyto_` | Dùng được cho |
|---|---|---|---|---|---|---|
| **all7** — đủ ở cả 7 cohort | **34** | 9 | 2 | 8 | 15 | **5.078 BN, 0 ô NA do panel** |
| **ge5** — ≥ 5 cohort | 33 | 24 | 1 | 3 | 5 | 4.000–5.000 BN tuỳ trường |
| **ge3** — ≥ 3 cohort | 88 | 79 | 9 | 0 | 0 | chủ yếu gene chỉ có ở UK-NCRI + 2 cohort WES |
| **cohort_specific** | 55 | 7 | 2 | 3 | 43 | chỉ dùng trong nội bộ một cohort |

**34 trường tầng all7:** `mut_ASXL1`, `mut_DNMT3A`, `mut_EZH2`, `mut_IDH1`, `mut_NPM1`, `mut_RUNX1`,
`mut_SRSF2`, `mut_TET2`, `mut_TP53` · `var_FLT3_ITD`, `var_FLT3_TKD` · `fus_RUNX1_RUNX1T1`,
`fus_CBFB_MYH11`, `fus_PML_RARA`, `fus_KMT2A_rearranged`, `fus_KMT2A_MLLT3`, `fus_DEK_NUP214`,
`fus_GATA2_MECOM`, `fus_BCR_ABL1` · 15 trường `cyto_` (complex, del_chr3/5/7/12/17/18, del_9q, del_20q,
gain_chr8/11/13/21/22, loss_chrY).

---

## 2. CÁC BƯỚC XỬ LÝ — số liệu và file output

| Bước | Nội dung |
|---|---|
| B0 | Thu thập nguồn |
| B1 | Đưa về 1 dòng = 1 bệnh nhân, **mẫu chẩn đoán** | 
| B2 | Chuẩn hoá thuật ngữ cột gene → trường codebook | 
| B3 | Mã hoá genotype 3 trạng thái + mask `assayed_` | 
| B4 | Chuẩn hoá outcome (OS, CR) |
| B5 | Khử trùng lặp bệnh nhân |
| B6 | Gộp 7 cohort + sinh ma trận 210 trường |
| B7 | Xử lý hiệu ứng cohort |

### B0 · Thu thập nguồn

| Cohort | File mutation | File outcome |
|---|---|---|---|
| UK-NCRI | `UK-NCRI.tsv` — 84 cột gene/variant + 70 cột cytogenetic, ma trận 0/1 | `os` + `os_status` (cùng file)|
| AMLSG | `AMSLG.tsv` — 81 cột, ma trận 0/1| `OS` + `Status`; `Post_Induction_Therapy_Remission` ở `AMLSG_Clinical_Anon.xlsx` |
| BeatAML 2.0 | `beataml_wes_wv1to4_mutation.csv` (11.721 variant / 871 mẫu WES) | `overallSurvival` + `vitalStatus`; `responseToInductionTx` |
| ClinSeq | `ClinSeq_Curated_mutations_AML.csv` (879 variant / 295 BN) + panel 25 gene | `follow_up_days` + `status`; `cr` | nội bộ |
| GSE146173 | 20 cột mutation trong series matrix | `overall_survival_(days)` + `life_status`; `resonse_to_induction_treatment` |
| TCGA-LAML | MAF `genome.wustl.edu_LAML...maf` (2.585 variant / 197 BN) + `SuppTable01.xlsx` | `days_to_death` / `days_to_last_followup` + `vital_status` |
| FIMM | `File_6_Binary_mutation_225S_57G.xlsx` | `os_time` + `os_event`; `best_disease_status` |

### B1 · Một dòng = một bệnh nhân, chỉ mẫu chẩn đoán

| Cohort | BN trong file nguồn | Mẫu chẩn đoán | 1 mẫu/BN | Có mutation | **Vào nhánh** |
|---|---|---|---|---|---|
| UK-NCRI | 2.113 | 2.113 | 2.113 | 2.113 | **2.113** |
| AMLSG | 1.540 | 1.540 | 1.540 | 1.540 | **1.540** |
| BeatAML | 805 *(942 dòng mức mẫu)* | **562** | 562 | 544 | **544** |
| ClinSeq | 315 | 315 | 315 | 315 | **315** |
| GSE146173 | 246 | 246 | 246 | 246 | **246** |
| TCGA-LAML | 200 | 200 | 197 | 197 | **186** |
| FIMM | 186 | **133** | 123 | 123 | **120** |
| **Tổng** | **5.405** | **5.109** | **5.096** | **5.078** | **5.064** |

Cột "vào nhánh" = thêm điều kiện có ≥ 1 outcome: TCGA-LAML mất 11 BN (không có `days_to_death` lẫn
`days_to_last_followup`), FIMM mất 3 BN. Bộ lọc thời điểm chẩn đoán:

- **BeatAML**: `diseaseStageAtSpecimenCollection == "Initial Diagnosis"` ∧ `isRelapse ≠ TRUE`
  (805 → 562 BN; loại 243 BN chỉ có mẫu Residual 202 / Relapse 63 / Remission 47 / Unknown).
  Sau đó `exomeSeq = y` → 544.
- **FIMM**: `Disease status == "Diagnosis"` trong `FIMM_sample_annot` (186 → 133 BN), rồi phải có cột
  trong `File_6` → 123.
- **TCGA-LAML**: toàn bộ là mẫu chẩn đoán (primary blood derived cancer).
- **UK-NCRI, AMLSG, ClinSeq, GSE146173**: mẫu lấy lúc vào thử nghiệm/chẩn đoán, không có mẫu tái phát.


### B2 · Chuẩn hoá thuật ngữ — tên cột gốc → trường codebook
**Viết codebook đích: mỗi trường một dòng, mỗi dòng có nguồn chuẩn.** Nguồn chọn theo khối:

| Khối | Nguồn |
|---|---|
| `mut_` / `fus_` / `gene_symbol` | **HGNC** (symbol hiện hành, previous/alias symbol) |
| `out_cr_*`, `clin_eln2017` | **ELN 2017** (Döhner 2017) — tiêu chí đáp ứng và phân tầng |
| `cyto_` | **ISCN** (chuỗi karyotype) |
| `clin_who2022`, `clin_aml_onset` | **WHO 2022 / ICC 2022** |
| `clin_wbc/hb/platelets/blast` | đơn vị **SI** |
| `out_mrd_*` | **ELN MRD 2021** |

→ Output: codebook dạng CSV + YAML (**243 trường**, core 53 / extended 190).
→ Kiểm: mọ

**mapping_dictionary** Ánh xạ các trường trong cohort vào codebook.

| Số đo | Giá trị |
|---|---|
| Cột gốc UK-NCRI đã map | 154 (84 gene/variant + 70 cytogenetic) → **153 trường đích** |
| Cột gốc AMLSG đã map | 81 (57 + 24) → **80 trường đích** |
| Cột chung hai cohort | 65 (riêng gene/variant **42/99**) |

### B3 · Mã hoá genotype 3 trạng thái + mask `assayed_`

Quy tắc theo loại xét nghiệm:

| Loại nguồn | Cohort | Cách sinh giá trị | Cách sinh mask |
|---|---|---|---|
| Ma trận 0/1 | UK-NCRI, AMLSG | dùng trực tiếp | `assayed = 1` cho đúng các cột cohort phát hành, **NaN** cho phần còn lại |
| WES/WGS (danh sách variant) | BeatAML, TCGA-LAML | pivot variant → wide; **loại biến thể silent/intron/UTR** (giữ missense, nonsense, frameshift, splice site, in-frame indel, start/stop lost) | `assayed = 1` cho **toàn bộ 119 gene** của universe (giải trình tự toàn hệ gene) |
| Panel đích danh | ClinSeq (25 gene), FIMM (57), GSE146173 (16) | pivot bảng curated / ma trận nhị phân | `assayed = 1` **chỉ cho gene trong panel**, NaN cho gene ngoài panel |

Tỉ lệ ô được xét nghiệm theo cohort: TCGA-LAML 76,6 % · BeatAML 74,8 % · UK-NCRI 73,3 % · FIMM 42,5 % ·
AMLSG 39,0 % · ClinSeq 30,5 % · GSE146173 23,0 %. Toàn ma trận **57,4 %**.

### B4 · Chuẩn hoá outcome

**OS** quy về **ngày**: UK-NCRI và AMLSG gốc là **năm** (× 365,25); AMLSG dùng cột `OS` của file clinical
vốn đã là ngày. Trung vị OS: AMLSG 774 ·
GSE146173 630 · UK-NCRI 573 · FIMM 415 · TCGA-LAML 365 · ClinSeq 327 · BeatAML 308 ngày.

**CR** mapping rồi mới nhị phân hoá (`out_cr_binary_strict`: CR/CRi = 1; RD và early_death = 0;
PR / not_evaluable / no_treatment / unknown = NA):

| Cohort | n có nhãn | % CR | Nguồn và quy tắc |
|---|---|---|---|
| UK-NCRI | 2.017 | **72,4 %** | `CR_stat`; 0 gộp cả refractory lẫn tử vong không đạt CR ⇒ chỉ có strict (CR-02) |
| AMLSG | 1.527 | **70,6 %** | `Post_Induction_Therapy_Remission` (CR 1.078 / RD 326 / ED 123 / PR 1 / NA 12) |
| BeatAML | 397 | 67,3 % | `responseToInductionTx`; `cumulativeChemo = n` → `no_treatment` |
| GSE146173 | 246 | 63,4 % | 6 trạng thái → CR/CRi/RD/early_death |
| ClinSeq | 308 | 57,8 % | cột `cr` yes/no |
| FIMM | 97 | 85,6 % | `best_disease_status` giải mã theo codebook `File_1.2` (1–4 → CR, 5 → CRi, 6 → PR, 7 → not_evaluable, 8/96 → RD, 94 → early_death, 97 → no_treatment) |
| TCGA-LAML | **0** | — | không có nguồn công khai (TCGA-01) |


### B5 · Khử trùng lặp bệnh nhân

| Cặp cohort | Kết quả |
|---|---|
| UK-NCRI ∩ AMLSG | **0 BN trùng** (giao 2.113 với 1.540 Sanger PD ID) |
| GSE106291 ∩ GSE146173 | 237/250 trùng — GSE106291 **không** đưa vào nhánh |


### B6 · Gộp và sinh ma trận

Mỗi cohort được đưa về cùng 210 cột; ô thuộc trường cohort không xét nghiệm để **NaN**; sau đó mask
`assayed_X = X.notna()`. Kết quả: **5.078 × (210 + 210)**.

### B7 · Hiệu ứng cohort — phương án đã thực hiện

Ba lớp xử lý:

1. **Cố định feature set theo độ phủ panel**  dùng tầng **all7** (34 trường, 0 ô NA do panel).
2. **`source_cohort` là biến phân tầng**, giữ nguyên trong `aml_core.csv`; `qc_split` được phân tầng theo
   cohort (80/20, seed 42) để train/test có cùng cơ cấu cohort.
3. **Đánh giá kết quả** (mục 4).

---

## 3. TỪNG COHORT GỐC — xử lý, mapping, làm sạch

| Cohort | BN nguồn → vào nhánh | OS | CR | Số trường được xét nghiệm | % ô được xét nghiệm |
|---|---|---|---|---|---|
| UK-NCRI | 2.113 → **2.113** | 2.113 | 2.017 | 154 | 73,3 % |
| AMLSG | 1.540 → **1.540** | 1.540 | 1.527 | 82 | 39,0 % |
| BeatAML | 805 → **544** | 544 | 397 | 161 | 74,8 % |
| ClinSeq | 315 → **315** | 315 | 308 | 65 | 30,5 % |
| GSE146173 | 246 → **246** | 246 | 246 | 50 | 23,0 % |
| TCGA-LAML | 200 → **186** | 186 | 0 | 161 | 76,6 % |
| FIMM | 186 → **120** | 120 | 97 | 92 | 42,5 % |

### 3.1. UK-NCRI (Tazi 2022)
Ma trận genotype đặc hoàn toàn (0 ô NA) ⇒ mọi BN vào nhánh. OS tính bằng năm → × 365,25. CR ở file supplement riêng, phủ 2.017/2.113. Không phát hành
chuỗi karyotype ISCN — 70 cột cytogenetic đã mã hoá sẵn là nguồn duy nhất, nên **`cyto_normal` = NA** cho toàn cohort.

### 3.2. AMLSG
`AMSLG.tsv` phân tách bằng khoảng trắng; khớp 1:1 PDID với file clinical (1.540/1.540). CR dạng text 4 mức
nên tách được đầy đủ `early_death` ⇒ có cả `out_cr_binary` lẫn `out_cr_binary_strict`. `cyto_normal` lấy
từ cột `NK` của file clinical.

### 3.3. BeatAML 2.0
File clinical ở mức mẫu (942 dòng → 805 subject). Lọc **chỉ mẫu chẩn đoán** còn 562 BN; thêm `exomeSeq = y`
còn **544**. `karyotype` dạng ISCN có ở 494/544 BN; `consensusAMLFusions` (từ RNA-seq) kết hợp với karyotype để sinh khối `fus_`.

### 3.4. ClinSeq
315/315 BN có panel (`ClinSeq_Status.panel` = FB/FBx); 295 BN có ≥ 1 variant curated, 20 BN còn lại là
**âm tính toàn panel** (assayed = 1, giá trị 0), không phải Unknown. Panel chỉ 25 gene ⇒ 30,5 % ô được
xét nghiệm — thấp thứ hai. Giới tính: mã 1 = Male, 2 = Female (khớp 315/315 với `AML_merged_labels.csv`).

### 3.5. GSE146173
Mutation nhúng trong series matrix⇒ giữ đủ 246 BN. Cytogenetics là **nhãn nhóm** (`CN-AML`, `t(8;21)`, `Complex
Cytogenetics`…) chứ không phải chuỗi ISCN ⇒ chỉ suy được một phần khối `cyto_`/`fus_`; 13 BN ghi
`Various` được coi là **đã xét nghiệm nhưng không thuộc nhóm nào** (gán giá trị 0).

### 3.6. TCGA-LAML
197/200 BN có trong MAF; thêm 11 BN không có bất kỳ mốc thời gian sống còn nào ⇒ **186** vào nhánh.
Không có CR (TCGA-01). Fusion lấy từ các cột curated của `SuppTable01` (PML-RARA, MYH11-CBFB,
RUNX1-RUNX1T1, BCR-ABL, MLL-partner, NUP98-NSD1) và `MLL-PTD` → `var_KMT2A_PTD`; karyotype ISCN từ cột
`Cytogenetics`.

### 3.7. FIMM 
186 BN → **133** có mẫu chẩn đoán → **123** có cột trong `File_6` → **120** có outcome (−35 %).
Nguyên nhân: nhiều BN chỉ được giải trình tự ở lần tái phát. `best_disease_status` giải mã theo codebook
`File_1.2`.

---

## 4. NGHIỆM THU


### 4.1. Tần suất đột biến

| Cohort | NPM1 | FLT3-ITD | DNMT3A | TP53 | RUNX1 | Complex karyotype |
|---|---|---|---|---|---|---|
| UK-NCRI | 32 % | 21 % | 28 % | 8 % | 15 % | 10 % |
| AMLSG | 28 % | 24 % | 23 % | 6 % | 9 % | 10 % |
| BeatAML | 24 % | 22 % | 20 % | 10 % | 11 % | 19 % |
| ClinSeq | 29 % | 24 % | 22 % | 10 % | 15 % | 14 % |
| GSE146173 | 33 % | 20 % | 32 % | 8 % | 20 % | 9 % |
| TCGA-LAML | 27 % | 20 % | 26 % | 8 % | 9 % | 14 % |
| FIMM | 24 % | 21 % | 17 % | 3 % | 7 % | 9 % |

Tất cả nằm trong khoảng cho AML người lớn (NPM1 ~30 %, FLT3-ITD ~20–25 %, DNMT3A ~25 %, TP53 ~8 %).

---

**Tool:** pandas 2.x, scikit-learn, `lifelines` 0.30.3 (Cox PH cho LOCO). Nguồn chuẩn: HGNC (symbol),
ELN 2017 (Döhner 2017) cho nhãn CR và phân tầng nguy cơ, ISCN cho karyotype.

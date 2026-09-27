# GSE公開Rコード更新手順（再選抜型ブートストラップ追加版）

## 1. 既存コードで変更しないもの

既存の `01_GSE_main_and_diagnostics.R` ～ `05_GSE_make_figures.R` は、そのまま残してください。特に `01` は既に検証した中核計算コードなので、再選抜型ブートストラップを入れるために書き換えません。

今回追加する正式な手続きレベル再選抜ブートストラップは、`06_GSE_procedure_level_reselection_bootstrap.R` として分離します。

`01` に既にある inner reselection は、分散成分推定によるEBV変動と方策形成の同一標本依存を調べる**機構診断**です。`06` は各内側標本で REML -> EBLUP -> top-k選抜をやり直して basic 型95%区間そのものを評価する**手続きレベルの再選抜ブートストラップ**です。両者は別解析として残します。

## 2. この更新パッケージで追加・置換するファイル

既存リポジトリのルートに、次を追加してください。

- `06_GSE_procedure_level_reselection_bootstrap.R`：新規追加
- `README.md`：この版で置換

既存リポジトリの `validation/` フォルダでは、次を置換・追加してください。

- `validation_parse_all.R`：置換
- `run_small_self_tests.R`：置換
- `validate_06_reselection_bootstrap.R`：新規追加
- `validate_all_production_outputs.R`：置換
- `VALIDATION_REPORT.md`：置換

## 3. 最初に静的・小規模検証を行う

Rを起動する前に、ファイル名と配置を確認してください。その後、リポジトリのルートで次を実行します。

```bash
Rscript validation/validation_parse_all.R
Rscript validation/run_small_self_tests.R
Rscript validation/validate_06_reselection_bootstrap.R smoke
```

3つともエラーなく終了してから本解析へ進みます。

## 4. 既存の bootstrap master RDS を用意する

`06` は、`01` の `run_bootstrap_generator_diagnostic_final()` が作成した

```text
Simulation_II_bootstrap_generator_ALL.rds
```

を入力に使います。既に論文の固定方策型ブートストラップを作成したときの master RDS があるなら、再生成する必要はありません。

READMEの標準配置では、例えば

```text
paper_output/bootstrap_diagnostic/Simulation_II_bootstrap_generator_ALL.rds
```

です。

## 5. 外側標本の完全再構成を確認する

これは強く推奨します。内側B=500を回す前に、保存したmaster RDSから外側標本を再構成し、既存のConditionalとTrue-VCの結果が一致することを確認します。

```bash
Rscript validation/validate_06_reselection_bootstrap.R outer \
  paper_output/bootstrap_diagnostic/Simulation_II_bootstrap_generator_ALL.rds
```

`Exact outer reconstruction validation: PASS` になることを確認してください。

この処理は9条件×1,000外側反復のREML再計算を含むため、small testより時間がかかりますが、内側B=500のフル再選抜ブートストラップよりは軽いです。

## 6. パイロットを実行する

新しいRセッションで次を実行します。

```r
source("01_GSE_main_and_diagnostics.R", encoding = "UTF-8")
source("06_GSE_procedure_level_reselection_bootstrap.R", encoding = "UTF-8")

boot1000 <- readRDS(
  "paper_output/bootstrap_diagnostic/Simulation_II_bootstrap_generator_ALL.rds"
)

rboot20 <- run_reselection_bootstrap_from_master(
  boot1000 = boot1000,
  outer_indices = 1:20,
  B_inner = 100,
  output_dir = "paper_output/reselection_bootstrap_pilot",
  checkpoint_every = 5,
  resume = FALSE
)

print_reselection_bootstrap_from_master(rboot20)
```

パイロットは動作確認用です。20外側反復/B=100の数値を、論文のS=1,000/B=500の数値と直接比較して一致を要求しないでください。

## 7. 論文用フル解析を実行する

```r
rboot1000 <- run_reselection_bootstrap_from_master(
  boot1000 = boot1000,
  B_inner = 500,
  output_dir = "paper_output/reselection_bootstrap",
  checkpoint_every = 10,
  resume = TRUE,
  tolerance = 1e-8,
  stop_on_full_validation_failure = TRUE
)

print_reselection_bootstrap_from_master(rboot1000)
print(rboot1000$validation_table)
stopifnot(all(rboot1000$validation_table$within_tolerance %in% TRUE))
```

フルS=1,000/B=500では、少なくとも次の再現チェックがすべてPASSであることが必要です。

- outer Conditionalが保存済み結果と一致
- outer True-VCが保存済み結果と一致
- 再計算した固定方策型ブートストラップが保存済み論文値と一致

1つでもFAILなら、新しい再選抜結果を論文値として採用せず、RDS、乱数順序、checkpoint、R/BLAS環境を確認してください。既存値に合わせるためにコードを書き換えないでください。

## 8. 補足表S10・S11を再生成する

フル解析が終わった後は、保存されたper-outer診断から再計算できます。内側ブートストラップを再実行する必要はありません。

```r
rbootm_write_manuscript_diagnostic_tables(
  rboot1000,
  output_dir = "paper_output/reselection_bootstrap"
)
```

または、Rセッションを閉じた後なら、master RDSを読み直して実行します。

```r
source("01_GSE_main_and_diagnostics.R", encoding = "UTF-8")
source("06_GSE_procedure_level_reselection_bootstrap.R", encoding = "UTF-8")

rboot1000 <- readRDS(
  "paper_output/reselection_bootstrap/reselection_bootstrap_master_only_ALL.rds"
)

rbootm_write_manuscript_diagnostic_tables(
  rboot1000,
  output_dir = "paper_output/reselection_bootstrap"
)
```

## 9. 最終結果を論文値と照合する

```bash
Rscript validation/validate_06_reselection_bootstrap.R validate \
  paper_output/reselection_bootstrap/reselection_bootstrap_master_only_ALL.rds
```

`Procedure-level reselection-bootstrap production validation: PASS` を確認してください。

このvalidatorは、埋め込みの完全再現チェックに加え、本文表9と補足表S10・S11の主要数値が最終論文値と整合するかを確認します。

## 10. 全解析をまとめて検証する場合

既存のvalidatorに `--reselection` を追加しました。

```bash
Rscript validation/validate_all_production_outputs.R \
  --main paper_output/main/Simulation_II/Simulation_II_method_performance.csv \
  --truevc paper_output/trueVC/Simulation_II_true_VC_oracle_method_performance.csv \
  --rlup RLUP_prior_sensitivity_S1000_RNG_FIXED/RLUP_prior_sensitivity_MASTER.rds \
  --matched paper_output/RAMN_matched \
  --pedigree paper_output/GSE_pedigree_structure_sensitivity_production_S1000_B500/pedigree_structure_sensitivity_FINAL_comparison.csv \
  --reselection paper_output/reselection_bootstrap/reselection_bootstrap_master_only_ALL.rds
```

## 11. GitHub / Zenodo公開前

最後に以下を行ってください。

1. `validation_parse_all.R` と `run_small_self_tests.R` がPASS
2. `validate_06_reselection_bootstrap.R validate` がPASS
3. 既存のproduction validatorもPASS
4. definitive runの `sessionInfo()` を保存
5. `README.md` に論文DOIを追記（採択後）
6. software licenseを追加
7. GitHubでversion tag/releaseを作成
8. そのreleaseをZenodoへarchive
9. Zenodo DOIを論文のData/code availabilityへ記載

## 補足：今回変更しないもの

再選抜型ブートストラップの追加を理由に、既存の `01` のvalidated numerical logic、RL-UPのRNG順序、same/cross、段階置換、選抜強度感度などを変更する必要はありません。今回の更新は、最終論文で新たに使う手続きレベル再選抜ブートストラップを公開可能な形で追加し、その再現・検証経路をREADMEとvalidationへ組み込むためのものです。

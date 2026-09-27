# Windows/RStudio: [8/8] pedigree design のパスエラー修正

`validation/run_small_self_tests.R` から `validate_04_pedigree_structure.R` を介して
`04_GSE_pedigree_structure_sensitivity.R` を入れ子で `source()` した場合、R の
`sys.frame(1)$ofile` が validation 側のファイルを指すことがあり、旧版は
`validation/01_GSE_main_and_diagnostics.R` を探して失敗する場合がありました。

修正版 `04_GSE_pedigree_structure_sensitivity.R` は、必要な 01 の関数が既に読み込まれて
いれば再 source しません。未読込の場合だけ、(1) 現在のスクリプト位置、(2) その親、
(3) 作業ディレクトリから `01_GSE_main_and_diagnostics.R` を探します。

科学計算ロジックは変更していません。実行経路・ファイル探索だけの修正です。

## 再確認

release root を作業ディレクトリにして、RStudio Console から:

```r
source("validation/run_small_self_tests.R")
```

最後に

```text
ALL SMALL GSE PUBLIC-RELEASE SELF-TESTS PASSED
```

と出れば small self-tests は完了です。

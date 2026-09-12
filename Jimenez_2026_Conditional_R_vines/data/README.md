# Kenneth R. French data

The financial and permutation-entropy applications use the daily
value-weighted returns of the ten industry portfolios from the Kenneth R.
French Data Library:

```text
https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip
```

The article uses observations from 2015-01-01 through 2024-12-31 and retains
5,032 complete daily observations. The downloaded archive used for the frozen
analysis has:

```text
SHA-256: b1a23edba335ddd37223c46d75a6aa07827170dd6fd04d9a27249602e1054e4f
Size: 933519 bytes
```

The public archive is not redistributed because its redistribution terms were
not independently verified. Download and verify it with:

```bash
./data/download_kenneth_french.sh
```

The script places verified copies at the two locations expected by the frozen
financial and ordinal scripts. Processing details are recorded in
`code/reproduction/financial/DATA_PROVENANCE.md` and
`code/reports/jmva_extension/permutation_entropy/final/DATA_PROVENANCE.md`.


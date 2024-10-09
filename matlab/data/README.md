# ct_dsge/data

This folder contains empirical data for the replication of the project.
- [dataset.xlsx](dataset.xlsx) contains U.S. on aggregate consumption, the average amount of hours worked and the interest rate for the period 1960:Q1 through
2019:Q4. All variables are obtained from the Federal Reserve Economic Data (FRED) database.
    -  **Aggregate consumption** is measured by monthly nominal personal consumption
      expenditures (`PCE`), deflated by the corresponding monthly price index (`PCEPI`).
      - Consumption at quarterly frequency is computed by aggregating monthly real expenditures over the quarter.
    - **Quarterly hours worked** correspond to the number of hours worked by wage and salary workers on nonfarm payrolls (`TOTLQ`).
    - All variables are transformed into per-capita values using the civilian non-institutional population aged 16 and over (`CNP16OV`) from the U.S. Bureau of Labor Statistics. With the exception of population, all variables are seasonally adjusted.
    - **Nominal 3-month Treasury bill rate** (`TB3MS`) from FRED, deflated using the `PCEPI`. Observations are end-of-quarter values. 

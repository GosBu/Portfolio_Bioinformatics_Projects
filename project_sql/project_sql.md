# pd4657-DATAB-projekt
## Wyniki zapytań SQL

### Zapytanie 2: Testy z więcej niż jednym pacjentem
**Zapytanie:**
```sql

SELECT tests.test_type, COUNT(patient_test.patient_id) AS patient_count
    FROM tests
    INNER JOIN patient_test ON tests.test_id = patient_test.test_id
    GROUP BY tests.test_type
    HAVING COUNT(patient_test.patient_id) > 30;
```

**Wynik**
```
+-----------+---------------+
| test_type | patient_count |
+-----------+---------------+
| SNP Array |            49 |
| NGS-panel |            43 |
+-----------+---------------+
```

### Zapytanie 3: Testy z więcej niż jednym pacjentem
**Zapytanie:**
```sql

SELECT tests.test_type, results.variant, COUNT(*) AS variant_count
    FROM results
    JOIN tests ON results.test_id = tests.test_id
    GROUP BY tests.test_type, results.variant
    ORDER BY tests.test_type, variant_count DESC;
```
**Wynik**
```
+-----------+----------+---------------+
| test_type | variant  | variant_count |
+-----------+----------+---------------+
| NGS-panel | rs67890  |            29 |
| NGS-panel | rs54321  |            29 |
| NGS-panel | rs112233 |            29 |
| NGS-panel | rs98765  |            22 |
| NGS-panel | rs12345  |            20 |
| SNP Array | rs12345  |            37 |
| SNP Array | rs54321  |            35 |
| SNP Array | rs67890  |            31 |
| SNP Array | rs98765  |            24 |
| SNP Array | rs112233 |            23 |
| WES       | rs67890  |            27 |
| WES       | rs54321  |            27 |
| WES       | rs98765  |            27 |
| WES       | rs112233 |            24 |
| WES       | rs12345  |            15 |
| WGS       | rs54321  |            28 |
| WGS       | rs67890  |            21 |
| WGS       | rs98765  |            19 |
| WGS       | rs12345  |            18 |
| WGS       | rs112233 |            15 |
+-----------+----------+---------------+
```

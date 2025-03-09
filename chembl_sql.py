import psycopg2
import pandas as pd

# Define connection parameters
conn_params = {
    "host": "mem2",
    "database": "chembl",
    "user": "chembl",
    "port": 5433
    # If your database requires a password, include it here, e.g., "password": "your_password"
}

# SQL query to be executed
query = """
SELECT DISTINCT 
    s.canonical_smiles,
    act.pchembl_value,
    cs.accession AS uniprot_id
FROM compound_structures s
JOIN molecule_dictionary m ON s.molregno = m.molregno
JOIN compound_records r ON m.molregno = r.molregno
JOIN activities act ON r.record_id = act.record_id
JOIN assays a ON act.assay_id = a.assay_id
JOIN target_dictionary t ON a.tid = t.tid
JOIN target_components tc ON t.tid = tc.tid
JOIN component_sequences cs ON tc.component_id = cs.component_id
WHERE act.standard_type IN ('Ki', 'IC50')
  AND act.pchembl_value BETWEEN 4 AND 9;
"""

try:
    # Establish a connection to the database
    with psycopg2.connect(**conn_params) as conn:
        # Use pandas to execute the query and load the results into a DataFrame
        df = pd.read_sql_query(query, conn)

    # Save the DataFrame to a CSV file without the index column
    output_csv = "query_results.csv"
    df.to_csv(output_csv, index=False)
    print(f"Query results have been successfully saved to '{output_csv}'.")

except Exception as e:
    print("An error occurred while executing the query or saving the CSV file:")
    print(e)
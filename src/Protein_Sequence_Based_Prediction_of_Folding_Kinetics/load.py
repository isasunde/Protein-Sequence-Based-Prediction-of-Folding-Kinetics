import pandas as pd
from mysql import connector
from mysql.connector import Error



def load_data(data, connection, table_name, columns):
    """
    Load a cleaned DataFrame into an existing MySQL table.

    Parameters
    ----------
    data : pandas.DataFrame
        Cleaned dataset.
    connection : mysql.connector.MySQLConnection
        Active MySQL connection.
    table_name : str
        Destination table name.
    columns : dict
        Mapping from source column names to SQL column names.

    Returns
    -------
    int
        Number of rows inserted.
    """
    if data.empty:
        print("No records to load.")
        return 0

    missing_columns = set(columns) - set(data.columns)
    if missing_columns:
        raise ValueError(
            f"Missing required DataFrame columns: {sorted(missing_columns)}"
        )

    # Convert pandas missing values (NaN, NaT, pd.NA) to Python None.
    data = data.astype(object).where(pd.notna(data), None)

    column_sql = ", ".join(f"`{col}`" for col in data.columns)
    placeholders = ", ".join(["%s"] * len(data.columns))

    query = (
        f"INSERT INTO `{table_name}` ({column_sql}) "
        f"VALUES ({placeholders})"
    )

    records = list(data.itertuples(index=False, name=None))
    cursor = connection.cursor()

    try:
        cursor.executemany(query, records)
        connection.commit()

        inserted = cursor.rowcount
        print(f"Successfully inserted {inserted} records into {table_name}.")
        return inserted

    except Error:
        connection.rollback()
        raise

    finally:
        cursor.close()

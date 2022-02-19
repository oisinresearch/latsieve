#include "sqlite3.h"
#include <stdio.h>
#include <stddef.h>
#include <cstdlib>
#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <iomanip> // setprecision
#include <gmpxx.h>
#include <cmath>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <cstring>	// memset
#include <omp.h>

using std::cout;
using std::endl;
using std::flush;
using std::string;
using std::ifstream;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::sort;
using std::to_string;
using std::hex;
using std::stringstream;
using std::abs;

int main(int argc, char** argv)
{
	if (argc != 4) {
		cout << endl << "Usage: ./specialqdb dbname q0max q1max" << endl << endl;
		return 0;
	}

	cout << "# ";
	for (int i = 0; i < argc; i++) cout << argv[i] << " ";
	cout << endl;

	string dbname = string(argv[1]);
	int64_t q0max = strtoll(argv[2], NULL, 10);
	int64_t q1max = strtoll(argv[3], NULL, 10);
	mpz_t qmpz0; mpz_init(qmpz0);
	mpz_t qmpz1; mpz_init(qmpz1);
	int64_t q0 = 1;
	int64_t q1 = 1;

    sqlite3* db = NULL;
    sqlite3_stmt* query1 = NULL;
    sqlite3_stmt* query2 = NULL;
    int ret = 0;
	string sql = "";
	// initialize engine
	if (SQLITE_OK != (ret = sqlite3_initialize())) {
		printf("Failed to initialize library: %d\n", ret);
		goto done;
	}

	// open connection to a DB
	if (SQLITE_OK != (ret = sqlite3_open_v2(dbname.c_str(), &db,
		SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL))) {
		printf("Failed to open conn: %d\n", ret);
		goto done;
	}

	// prepare the statement
	sql = "CREATE TABLE specialq0 (q INT NOT NULL, count INT, PRIMARY KEY (q));";
	if (SQLITE_OK != (ret = sqlite3_prepare_v2(db, sql.c_str(), -1, &query1, NULL))) {
		printf("Failed to prepare create table: %d, %s\n", ret, sqlite3_errmsg(db));
		goto done;
	}

	// step to 1st row of data
	if (SQLITE_DONE != (ret = sqlite3_step(query1))) {
		printf("Failed to create table: %d, %s\n", ret, sqlite3_errmsg(db));
		goto done;
	}

	// prepare the statement
	sql = "CREATE TABLE specialq1 (q INT NOT NULL, count INT, PRIMARY KEY (q));";
	if (SQLITE_OK != (ret = sqlite3_prepare_v2(db, sql.c_str(), -1, &query2, NULL)))
	{
		printf("Failed to prepare create table: %d, %s\n", ret, sqlite3_errmsg(db));
		goto done;
	}

	// step to 1st row of data
	if (SQLITE_DONE != (ret = sqlite3_step(query2))) {
		printf("Failed to create table: %d, %s\n", ret, sqlite3_errmsg(db));
		goto done;
	}

	// prepare the statement
	sql = "INSERT INTO specialq0 (q, count) VALUES (?, NULL);";
	if (SQLITE_OK != (ret = sqlite3_prepare_v2(db, sql.c_str(), -1, &query1, NULL))) {
		printf("Failed to prepare insert: %d, %s\n", ret, sqlite3_errmsg(db));
		goto done;
	}

	// prepare the statement
	sql = "INSERT INTO specialq1 (q, count) VALUES (?, NULL);";
	if (SQLITE_OK != (ret = sqlite3_prepare_v2(db, sql.c_str(), -1, &query2, NULL))) {
		printf("Failed to prepare insert: %d, %s\n", ret, sqlite3_errmsg(db));
		goto done;
	}

	// store in memory
	sqlite3_exec(db, "PRAGMA journal_mode = MEMORY", NULL, NULL, NULL);

	while (q0 < q0max || q1 < q1max) {
		mpz_set_ui(qmpz0, q0);
		mpz_nextprime(qmpz0, qmpz0);
		q0 = mpz_get_ui(qmpz0);
		if (q0 < q0max) {
			sqlite3_bind_int(query1, 1, q0);
			// execute query
			if (SQLITE_DONE != (ret = sqlite3_step(query1))) {
				printf("Failed to insert: %d, %s\n", ret, sqlite3_errmsg(db));
				goto done;
			}
			sqlite3_clear_bindings(query1);
			sqlite3_reset(query1);
		}
		mpz_set_ui(qmpz1, q1);
		mpz_nextprime(qmpz1, qmpz1);
		q1 = mpz_get_ui(qmpz1);
		if (q1 < q1max) {
			sqlite3_bind_int(query2, 1, q1);
			// execute query
			if (SQLITE_DONE != (ret = sqlite3_step(query2))) {
				printf("Failed to insert: %d, %s\n", ret, sqlite3_errmsg(db));
				goto done;
			}
			sqlite3_clear_bindings(query2);
			sqlite3_reset(query2);
		}
	}

	done:

    // cleanup
    if (NULL != query1) sqlite3_finalize(query1);
    if (NULL != query2) sqlite3_finalize(query2);
    if (NULL != db) sqlite3_close(db);
    sqlite3_shutdown();
	mpz_clear(qmpz0);
	mpz_clear(qmpz1);

    return 0;
}


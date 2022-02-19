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
#include <chrono>
#include <thread>

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
using std::this_thread::sleep_for;
using std::chrono::milliseconds;

#define BUFFER_SIZE 256

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
	mpz_t pmpz0; mpz_init(pmpz0);
	int64_t q0 = 65536;
	int64_t q1 = 65536;

    sqlite3* db = NULL;
    sqlite3_stmt* query1 = NULL;
    sqlite3_stmt* query2 = NULL;
	char * errmsg = 0;
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

	// wait up to 100ms to try transactions
	sqlite3_busy_timeout(db, 1000);

	// store in memory
	sqlite3_exec(db, "PRAGMA journal_mode = MEMORY", NULL, NULL, &errmsg);

	sql = "UPDATE specialq0 SET count=0 WHERE q=(SELECT MAX(q) FROM specialq0 "
		  "WHERE count IS NULL) RETURNING q;";
	sqlite3_prepare_v2(db, sql.c_str(), BUFFER_SIZE, &query1, NULL);

	sql = "UPDATE specialq0 SET count = CASE WHEN count IS NULL THEN 1 "
		  "ELSE count + 1 END WHERE q = ?;";
	sqlite3_prepare_v2(db, sql.c_str(), BUFFER_SIZE, &query2, NULL);

	srand(time(NULL));

	while (q0 > 10000) {
		// prepare the statement
		//sql = "begin exclusive transaction; update specialq0 set count=0 where "
		//	  "q=(select max(q) from specialq0 where count is NULL) returning q; "
		//	  "commit transaction;";
		// execute query
		if (SQLITE_ROW != (ret = sqlite3_step(query1))) {
			printf("Failed to SELECT MAX(q): %d, %s\n", ret, sqlite3_errmsg(db));
			goto done;
		}
		string qstr = string(reinterpret_cast<const char*>(sqlite3_column_text(query1, 0)));
		cout << qstr << endl;
		q0 = stoll(qstr);
		if (SQLITE_DONE != (ret = sqlite3_step(query1))) {
			printf("SELECT MAX(q) not done: %d, %s\n", ret, sqlite3_errmsg(db));
			goto done;
		}
		// now write 10 random primes to specialq0 table (really sieve base primes)
		sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errmsg);
		for (int i = 0; i < 10; i++) {
			int p0 = rand() % (q0max - 1000);
			mpz_set_ui(pmpz0, p0);
			mpz_nextprime(pmpz0, pmpz0);
			p0 = mpz_get_ui(pmpz0);
			sqlite3_bind_int(query2, 1, p0);
			//cout << sqlite3_expanded_sql(query2) << endl;
			if (SQLITE_DONE != (ret = sqlite3_step(query2))) {
				printf("Failed to UPDATE: %d, %s\n", ret, sqlite3_errmsg(db));
				goto done;
			}
			sqlite3_clear_bindings(query2);
			sqlite3_reset(query2);
		}
		if (SQLITE_OK != (ret = sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &errmsg))) {
			printf("Failed to COMMIT: %d, %s\n", ret, sqlite3_errmsg(db));
			goto done;
		}
		// simulate delay to give other threads a chance
		//sleep_for(milliseconds(rand()%20));
	}

	done:

    // cleanup
    if (NULL != query1) sqlite3_finalize(query1);
    if (NULL != query2) sqlite3_finalize(query2);
    if (NULL != db) sqlite3_close(db);
    sqlite3_shutdown();
	mpz_clear(pmpz0);
	mpz_clear(qmpz0);
	mpz_clear(qmpz1);

    return 0;
}


#ifndef _MCPELE_SQLITEDB_H
#define _MCPELE_SQLITEDB_H

#include <stdio.h>
#include <sqlite3.h>

namespace mcpele{

/*
 * SQL C interface crash course, for more details refer to  https://www.sqlite.org/cintro.html
 *
 * To run an SQL statement, the application follows these steps:
 *    1.  Create a prepared statement using sqlite3_prepare().
 *    2.  Evaluate the prepared statement by calling sqlite3_step() one or more times.
 *    3.  For queries, extract results by calling sqlite3_column() in between two calls to sqlite3_step().
 *    4.  Destroy the prepared statement using sqlite3_finalize().
 *
 * Convenience wrappers:
 * The sqlite3_exec() interface is a convenience wrapper that carries out all four of the above steps with
 * a single function call. A callback function passed into sqlite3_exec() is used to process each row of
 * the result set. This is also known as the One-Step Query Execution Interface. We try using this simplified
 * interface rather than constructing the 4 steps procedure outlined above. Details of sqlite3_exec:
 *
 * int sqlite3_exec(
 *   sqlite3*,                                  !An open database!
 *   const char *sql,                           !SQL to be evaluated!
 *   int (*callback)(void*,int,char**,char**),  !Callback function!
 *   void *,                                    !1st argument to callback!
 *   char **errmsg                              !Error msg written here!
 *   );
 *
 * If the callback function of the 3rd argument to sqlite3_exec() is not NULL, then it is invoked for each
 * result row coming out of the evaluated SQL statements. The 4th argument to sqlite3_exec() is relayed
 * through to the 1st argument of each callback invocation. If the callback pointer to sqlite3_exec() is NULL,
 * then no callback is ever invoked and result rows are ignored.
 * To avoid memory leaks, the application should invoke sqlite3_free() on error message strings returned through
 * the 5th parameter of of sqlite3_exec() after the error message string is no longer needed.
 *
 *  Callback function:
 *  If an sqlite3_exec() callback returns non-zero, the sqlite3_exec() routine returns SQLITE_ABORT without invoking
 *  the callback again and without running any subsequent SQL statements. The 2nd argument to the sqlite3_exec()
 *  callback function is the number of columns in the result. The 3rd argument to the sqlite3_exec() callback is an
 *  array of pointers to strings obtained as if from sqlite3_column_text(), one for each column. If an element of a
 *  result row is NULL then the corresponding string pointer for the sqlite3_exec() callback is a NULL pointer. The
 *  4th argument to the sqlite3_exec() callback is an array of pointers to strings where each entry represents the
 *  name of corresponding result column as obtained from sqlite3_column_name().
 *
 */

class CheckpointSqlite3{
private:
    sqlite3 *m_db; //database connection handle (pointer to sqlite3 structure)
    bool m_open;
    std::stringstream m_output;
public:
    CheckpointSqlite3();
    ~CheckpointSqlite3(){
        this->close();
    }

    /*sqlite_callback is a call back for which data is the 1st argument and
     * errmsg will be return to capture any error raised by the routine.
     *
     * This function will actually print a particular result, it will not return
     * a string or value that one can use, as is */
    static int callback(void *NotUsed, int argc, char **argv, char **azColName){
       for(int i=0; i<argc; i++){
           std::printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
       }
       std::printf("\n");
       return 0;
    }

    /*this is a wrapper for the constructor/destructor of the database connection handle
     * Each open SQLite database is represented by a pointer to an instance of
     * the opaque structure named "sqlite3". It is useful to think of an sqlite3
     * pointer as an object. sqlite3_open() amd sqlite3_close() are its constructor
     * and destructor respectively
     */

    //open sql database fname, if it doesn't exist it creates one
    void open(const char* fname)
    {
        int rc = sqlite3_open(fname, &m_db);
        if( rc ){
            std::fprintf(std::stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
            sqlite3_close(m_db);
            std::exit(0);
        }
        else{
            std::fprintf(std::stderr, "Opened database successfully\n");
        }
    }

    //close sql database
    void close(){
        sqlite3_close(m_db);
    }

    //create sql table
    //example of sql command to create table
    /* Create SQL statement
     *        sql_command = \
     *        "CREATE TABLE COMPANY("  \
     *        "ID INT PRIMARY KEY     NOT NULL," \
     *        "NAME           TEXT    NOT NULL," \
     *        "AGE            INT     NOT NULL," \
     *        "ADDRESS        CHAR(50)," \
     *        "SALARY         REAL );";
     */

    class Column{
        std::string _name; //eg column1
        std::string _type; //eg INT, TEXT, CHAR(50), KEY
        std::string _flags; //eg NOT NULL

        Column(const std::string name, const std::string type, const std::string flags):
        _name(name),_type(type),_flags(flags){}
        ~Column(){}

        //concatenate name type and flags to create a table column
        std::string get_cmd(){
            std::string cmd = _name + " " + _type + " " + _flags;
            return cmd;
        }
    };

    void create_table(std::string table_name, std::vector<Column> column_list){
        char *zErrMsg = 0;
        std::string sql_command = "CREATE TABLE " + table_name +"(";
        for(auto& column : column_list){
            sql_command += column.get_cmd() + ",";
        }
        sql_command.pop_back(); //remove last comma
        sql_command += ");";    //closing bracket
        int rc = sqlite3_exec(m_db, sql_command, callback, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            std::fprintf(std::stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
        else{
            std::fprintf(stdout, "Table created successfully\n"); //a command should be for more than creating a directory
        }
    }

    //example of how to create records in a table
    /* Create SQL statement
     *        sql = "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY) "  \
     *        "VALUES (1, 'Paul', 32, 'California', 20000.00 ); " \
     *        "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY) "  \
     *        "VALUES (2, 'Allen', 25, 'Texas', 15000.00 ); "     \
     *        "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY)" \
     *        "VALUES (3, 'Teddy', 23, 'Norway', 20000.00 );" \
     *        "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY)" \
     *        "VALUES (4, 'Mark', 25, 'Rich-Mond ', 65000.00 );";
    */

};

}
#endif




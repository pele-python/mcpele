#ifndef _MCPELE_SQLITEDB_H
#define _MCPELE_SQLITEDB_H

#include <stdio.h>
#include <sqlite3.h> //requires sudo apt-get install libsqlite3-dev
#include <iostream>
#include <sstream>
#include <cstdint>
#include <mcpele/serialize_vector.h>

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
    std::string m_db_name, m_callback_output;
    bool m_open;
public:

    CheckpointSqlite3(const std::string db_name):
        m_db(),
        m_db_name(db_name),
        m_callback_output(),
        m_open(false)
    {};

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

    /*assumes only one table and one column
    https://stackoverflow.com/questions/12984248/sqlite-in-c-c-sqlite3-exec-parameter-set-in-callback-function-is-pointing-to
    https://stackoverflow.com/questions/18839799/retrieve-sqlite-table-data-in-c
    //assumes a single value is being selected
    */
    static int callback_select(void *data, int num_fields, char **fields, char **azColName){
        std::string* result_str = static_cast<std::string*>(data);
        try {
            result_str->append(fields[0],num_fields); //DEBUG NOT SURE ABOUT THIS (AFRAID THIS IS GOING TO APPEND A SINGLE CHARACTER)
         }
         catch (...) {
           // abort select on failure, don't let exception propogate through sqlite3 call-stack
           return 1;
         }
        return 0;
    }

    /*this is a wrapper for the constructor/destructor of the database connection handle
     * Each open SQLite database is represented by a pointer to an instance of
     * the opaque structure named "sqlite3". It is useful to think of an sqlite3
     * pointer as an object. sqlite3_open() amd sqlite3_close() are its constructor
     * and destructor respectively
     */

    //open sql database fname, if it doesn't exist it creates one
    void open()
    {
        int rc = sqlite3_open(m_db_name.c_str(), &m_db);
        if( rc ){
            std::fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(m_db));
            sqlite3_close(m_db);
            std::exit(0);
        }
        else{
            std::fprintf(stderr, "Opened database successfully\n");
            m_open = true;
        }
    }

    //close sql database
    void close(){
        sqlite3_close(m_db);
    }

    void drop_table(std::string table_name){
        char *zErrMsg = 0;
        std::string sql_command = "DROP TABLE IF EXISTS " + table_name +";";
        int rc = sqlite3_exec(m_db, sql_command.c_str(), callback, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            std::fprintf(stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
        else{
            std::fprintf(stdout, "Table dropped successfully\n"); //a command should be for more than creating a directory
        }
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
    private:
        std::string _name; //eg column1
        std::string _type; //eg INT, TEXT, CHAR(50), KEY ...
        std::string _flags; //eg NOT NULL, UNIQUE, PRIMARY KEY, FOREIGN KEY, CHECK, DEFAULT

        Column(const std::string name, const std::string type, const std::string flags=""):
        _name(name),_type(type),_flags(flags){}
        ~Column(){}
    public:
        //concatenate name type and flags to create a table column
        std::string get_cmd(){
            std::string cmd = _name + " " + _type + " " + _flags;
            return cmd;
        }
    };

    //primary key is set by default to be ID (integer) by default
    //creates a table if it does not exists
    void create_table(std::string table_name, std::vector<Column> column_list, std::string primary_key="ID INT PRIMARY KEY NOT NULL"){
        char *zErrMsg = 0;
        std::string sql_command = "CREATE TABLE IF NOT EXISTS " + table_name +"(";
        sql_command += primary_key + ",";
        for(auto& column : column_list){
            sql_command += column.get_cmd() + ",";
        }
        sql_command.erase(sql_command.end() - 1); //remove last comma
        sql_command += ");";    //closing bracket
        int rc = sqlite3_exec(m_db, sql_command.c_str(), callback, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            std::fprintf(stderr, "SQL error: %s\n", zErrMsg);
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

    template<typename T>
    void insert(std::string table_name, std::string column_name, T value){
        static_assert(std::is_arithmetic<T>::value,"type is not arithmetic"); //checks that the type is integer or float
        char *zErrMsg = 0;
        std::string sql_command = "INSERT INTO " + table_name + " ";
        sql_command += "(" + column_name + ")";
        sql_command += "VALUES";
        sql_command += "(" + std::to_string(value) + ");";

        int rc = sqlite3_exec(m_db, sql_command.c_str(), callback, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            std::fprintf(stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
        else{
            std::fprintf(stdout, "Records created successfully\n");
        }
    }

    template<typename T>
    void insert_vector(std::string table_name, std::string column_name, std::vector<T> vec){
        char *zErrMsg = 0;
        std::string sql_command = "INSERT INTO " + table_name + " ";
        sql_command += "(" + column_name + ")";
        sql_command += "VALUES";
        std::string ser_vec = mcpele::serialize_vector<T>(vec);
        sql_command += "(" + ser_vec + ");";

        int rc = sqlite3_exec(m_db, sql_command.c_str(), callback, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            std::fprintf(stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
        else{
            std::fprintf(stdout, "Records created successfully\n");
        }
    }

    //select arithmetic value (there's a template specialisation for string)
    template<typename T>
    T select_value(std::string table_name, std::string column_name){
        static_assert(std::is_arithmetic<T>::value,"type is not arithmetic"); //checks that the type is integer or float
        char *zErrMsg = 0;
        std::string data;
        std::string sql_command = "SELECT " + column_name + " from " + table_name;

        int rc = sqlite3_exec(m_db, sql_command.c_str(), callback_select, &data, &zErrMsg);
        if( rc != SQLITE_OK ){
            std::fprintf(stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
        else{
            std::fprintf(stdout, "Operation done successfully\n");
        }

        return data;
    }

};
    //these should go in a separate c++ file
    //template specialisation for strings, these need to be outside the scope of the class
    template<>
    void CheckpointSqlite3::insert<std::string>(std::string table_name, std::string column_name, std::string value){
        char *zErrMsg = 0;
        std::string sql_command = "INSERT INTO " + table_name + " ";
        sql_command += "(" + column_name + ")";
        sql_command += "VALUES";
        sql_command += "(" + value + ");";

        int rc = sqlite3_exec(m_db, sql_command.c_str(), callback, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            std::fprintf(stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
        else{
            std::fprintf(stdout, "Records created successfully\n");
        }
    }

    //select arithmetic value (there's a template specialisation for string)
    template<>
    std::string CheckpointSqlite3::select_value<std::string>(std::string table_name, std::string column_name){
        char *zErrMsg = 0;
        std::string sql_command = "SELECT " + column_name + " from " + table_name;

        int rc = sqlite3_exec(m_db, sql_command.c_str(), callback_select, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            std::fprintf(stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
        else{
            std::fprintf(stdout, "Operation done successfully\n");
        }

        return m_callback_output;
    }

    //require specialisations for

}
#endif




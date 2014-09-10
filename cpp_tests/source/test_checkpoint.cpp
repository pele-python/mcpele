#include <iostream>
#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#include <cstdint>

#include "mcpele/sqlitedb.h"

TEST(CheckpointSqlite3, OpenClose){
    mcpele::CheckpointSqlite3 db("test.sqlite");
    db.open();
    db.close();
}

TEST(CheckpointSqlite3, Serialize_Deserialize_DBL){
    std::vector<double> vec;
    vec.push_back(0.);
    vec.push_back(1.);
    vec.push_back(2.);
    std::vector<double> cpvec(vec);
    std::string ser_vec = mcpele::serialize_vector<double>(vec);
    for(size_t i=0;i<3;++i){
        EXPECT_EQ(vec[i],cpvec[i]);
    }
    std::vector<double> newvec = mcpele::deserialize_vector<double>(ser_vec);
    EXPECT_EQ(newvec.size(),3);
    for(size_t i=0;i<3;++i){
        EXPECT_EQ(newvec[i],cpvec[i]);
    }
}

TEST(CheckpointSqlite3, Serialize_Deserialize_INT32_T){
    std::vector<int32_t> vec;
    vec.push_back(0);
    vec.push_back(1);
    vec.push_back(2);
    std::vector<int32_t> cpvec(vec);
    std::string ser_vec = mcpele::serialize_vector<int32_t>(vec);
    for(size_t i=0;i<3;++i){
        EXPECT_EQ(vec[i],cpvec[i]);
    }
    std::vector<int32_t> newvec = mcpele::deserialize_vector<int32_t>(ser_vec);
    EXPECT_EQ(newvec.size(),3);
    for(size_t i=0;i<3;++i){
        EXPECT_EQ(newvec[i],cpvec[i]);
    }
}

TEST(CheckpointSqlite3, Serialize_Deserialize_UINT32_T){
    std::vector<uint32_t> vec;
    vec.push_back(0);
    vec.push_back(1);
    vec.push_back(2);
    std::vector<uint32_t> cpvec(vec);
    std::string ser_vec = mcpele::serialize_vector<uint32_t>(vec);
    for(size_t i=0;i<3;++i){
        EXPECT_EQ(vec[i],cpvec[i]);
    }
    std::vector<uint32_t> newvec = mcpele::deserialize_vector<uint32_t>(ser_vec);
    EXPECT_EQ(newvec.size(),3);
    for(size_t i=0;i<3;++i){
        EXPECT_EQ(newvec[i],cpvec[i]);
    }
}


/*TEST(CheckpointSqlite3, OpenClose){
    mcpele::CheckpointSqlite3 db("test.sqlite");
    db.open();
    db.close();
}*/


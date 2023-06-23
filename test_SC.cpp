// test code for SC
// Psudo code now

#include "SimplicialComplex.hpp"
#include <catch2/catch.hpp>

TEST_CASE("lnk-case1", "[SC][lkn]")
{
    // auto V = {
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0}
    // }; // 4 vertices

    auto F = {
        {0,3,2},
        {0,1,3},
        {1,2,3}
    }; // 3 Faces

    // dump it to (Tri)Mesh
    Mesh m(F);

    // get the tuple point to V(0), E(01), F(013)
    long hash = 0;
    Tuple t(0, 2, 1, hash);

    SimplicialComplex lnk_0 = lnk(Simplex(t, 0), m);
    SimplicialComplex lnk_1 = lnk(Simplex(t.sw(0, m), 0), m);
    SimplicialComplex lhs = SC_intersect(lnk_0, lnk_1, m);
    SimplicialComplex lnk_01 = lnk(Simplex(t, 1), m);
    SimplicialComplex lnk_10 = lnk(Simplex(t.sw(0,m), 1), m);
    

    REQUIRE(lnk_0.get_size() == 5);
    REQUIRE(lnk_1.get_size() == 5);
    REQUIRE(lnk_01.get_size() == 1);
    REQUIRE(lhs.get_size() == 3);

    REQUIRE(lnk_01 == lnk_10);

    REQUIRE(link_cond(t, m) == false);
}


TEST_CASE("lnk-case2", "[SC][lnk]")
{
    // auto V = {
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0},
    //     {0,0,0}
    // }; // 6 vertices

    auto F = {
        {0,3,1},
        {0,1,2},
        {0,2,4},
        {2,1,5}
    }; // 4 Faces

    // dump it to (Tri)Mesh
    Mesh m(F);

    // get the tuple point to V(0), E(01), F(012)
    long hash = 0;
    Tuple t(0, 2, 1, hash);

    SimplicialComplex lnk_0 = lnk(Simplex(t, 0), m);
    SimplicialComplex lnk_1 = lnk(Simplex(t.sw(0, m), 0), m);
    SimplicialComplex lhs = SC_intersect(lnk_0, lnk_1, m);
    SimplicialComplex lnk_01 = lnk(Simplex(t, 1), m);
    SimplicialComplex lnk_10 = lnk(Simplex(t.sw(0,m), 1), m);
    

    REQUIRE(lnk_0.get_size() == 7);
    REQUIRE(lnk_1.get_size() == 7);
    REQUIRE(lnk_01.get_size() == 2);

    REQUIRE(lhs == lnk_01);
    REQUIRE(lnk_01 == lnk_10);

    REQUIRE(link_cond(t, m) == true);
}
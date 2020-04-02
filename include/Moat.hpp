#pragma once
#include <vector>
#include "../include/BaseSieve.hpp"
#include "../include/BlockSieve.hpp"
using namespace std;


// Gather data for Gaussian moat problem.
class OctantMoat {
private:
    double jumpSize;
    uint64_t normBound;
    bool verbose;
    vector<vector<bool>> sieveArray;
    vector<gint> nearestNeighbors, currentComponent;
    vector<vector<gint>> allComponents;

public:
    explicit OctantMoat(double, uint64_t = 0, bool = true);
    void setNearestNeighbors();
    void exploreComponent(int32_t, int32_t);
    uint32_t getComponentSize();
    gint getComponentMaxElement();
    void exploreAllComponents();
    vector<gint> getCurrentComponent();
    void printCurrentComponent();
    vector<vector<gint>> getAllComponents();
};


// Derived from BlockSieve
class VerticalMoat : public BlockSieve {
private:
    // Static variables.
    static bool verbose;
    static double jumpSize;
    static uint32_t realPart;
    static uint32_t blockSize, dx, dy;
    static uint64_t sievingPrimesNormBound;
    static vector<gint> sievingPrimes, nearestNeighbors;

    // Instance variables.
    uint32_t x, y;
    int32_t upperWallYPunch;
    uint64_t countVisited;
    int32_t farthestRight;

public:
    // Call this static setter method before any instances of this class are created.
    static void setStatics(int32_t, double, bool = true);
    static void findVerticalMoat();

    VerticalMoat(uint32_t, uint32_t);
    void callSieve();
    bool exploreAtGint(int32_t, int32_t, bool = false);
    bool exploreLeftWall();
    void exploreUpperWall();
    pair<int32_t, int32_t> getNextBlock();
};




// Also derived from BlockSieve
class SegmentedMoat : public BlockSieve {
private:
    // Static variables.
    static bool verbose;
    static double jumpSize;
    static uint32_t previousdy;
    static uint64_t blockSize;
    static uint64_t sievingPrimesNormBound;
    static vector<gint> sievingPrimes, nearestNeighbors;

    // The index of the outer vector determines which component number the inner
    // vector corresponds with.
    static vector<vector<gint>> leftBoundary;

    // Holding counts of component sizes. Individual components are indexed by
    // unsigned longs, the index of this vector.
    static vector<uint64_t> componentSizes;

    // Instance variables.
    uint32_t x, dx, dy;
    vector<vector<gint>> rightBoundary;
    // Holds status of components. If component has not propagated to right
    // boundary or merged with another component, it can be forgotten.
    vector<bool> hasComponentPropagated;
    vector<vector<uint32_t>> leftComponentLookUp;

public:
    static void setStatics(double, bool = true);
    static void setSievingPrimes();
    static uint64_t getCountMainComponent();

    SegmentedMoat(uint32_t, uint32_t, uint32_t);
    void callSieve();
    void exploreComponent(uint32_t, bool = true);
    void exploreLeftBoundary();
    void exploreRightBoundary();
    void runSegment();
    bool hasMainComponentPropagated();
};
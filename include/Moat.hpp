#pragma once
#include <vector>
#include "../include/BaseSieve.hpp"
using namespace std;


// Gather data for Gaussian moat problem.
class OctantMoat {
private:
    uint64_t normBound;
    double jumpSize;
    bool verbose;
    vector<vector<bool>> sieveArray;
    vector<gint> nearestNeighbors, currentComponent;
    vector<vector<gint>> allComponents;

public:
    OctantMoat(uint64_t, double, bool = true);
    void setNearestNeighbors();
    void exploreComponent(int32_t, int32_t);
    uint32_t getComponentSize();
    gint getComponentMaxElement();
    void exploreAllComponents();
    vector<gint> getCurrentComponent();
    void printCurrentComponent();
    vector<gint> getUnexplored();
    vector<vector<gint>> getAllComponents();
};


class VerticalMoat {
private:
    uint32_t realPart;
    double jumpSize;
    bool verbose;
    uint64_t normBound;
    int32_t x;
    int32_t y;
    int32_t dx;
    int32_t dy;
    // These two vectors could be made static, but they do also depend on realPart and jumpSize.
    vector<gint> sievingPrimes, nearestNeighbors;
    vector<gint> currentComponent;
    vector<vector<bool>> sieveArray;

public:
    VerticalMoat(uint32_t, double, bool = true);
    void setNearestNeighbors();
    void setSieveArray();
    void printSieveArray();
    int exploreAtGint(int32_t, int32_t);
    void exploreLeftWall();
    void exploreUpperWall();


};
/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);

void computeAcceleration(struct world* jello, struct point a[8][8][8]);
point computeSpringForce(point pA, point pB, double kHook, double restLength);
point computeDampingForce(point pA, point pB, point vA, point vB, double kDamp);

// Spring
void calculateSpring(struct world* jello, struct point a[8][8][8]
    , int i, int j, int k, int nextI, int nextJ, int nextK, double restLength);
void calculateShearSpring(struct world* jello, struct point a[8][8][8], int i, int j, int k);

point interpolateForce(struct world* jello, int i, int j, int k);

// Collision 
bool isCollideWithBox(point p, double boxWidth);
bool isCollideWithPlane(point p, double a, double b, double c, double d);

point computePenaltyForce(struct world* jello, point p, point velocity);
point computePlanePenaltyForce(struct world* jello, point p, point velocity);


// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);

#endif


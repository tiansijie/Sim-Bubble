#ifndef BUBBLE_SIM_H_INCLUDED
#define BUBBLE_SIM_H_INCLUDED

//#include "openGL_headers.h"
#include "math_headers.h"
#include "glm.hpp"
#include <vector>
#include "bubble.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <cassert>
#include <vector>
#include <iterator>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>  Traits;

typedef Traits::RT                                          Weight;
typedef Traits::Bare_point                                  Point;
typedef Traits::Weighted_point                              Weighted_point;

typedef CGAL::Regular_triangulation_3<Traits>               Rt;

typedef Rt::Vertex_iterator                                 Vertex_iterator;
typedef Rt::Vertex_handle                                   Vertex_handle;


class Bubbles_sim{
public:
	Bubbles_sim();
	~Bubbles_sim();	

	//void draw(const VBO& vbos);
	void update(double dt);
	
	//basic operation
	void insertBubble(glm::vec3 center, glm::vec3 velocity, float radius);
	Bubble getBubble(int i);
	glm::vec3 getPosition(int i);
	glm::vec3 getVelocity(int i);
	float getRadius(int i);
	void setBubble(int i);
	int size();
protected:
	// Topology methods
	void regularTriangulation_3();
	//forces function
	void StrongBubbleForce();
	void WeakBubbleForce();
	void DragForce();
	void BuoyantForce();
	void LiquidAdhesionForce();
	void SolidAdhesionForce();
	void SolidAttractionForce();
	void DampingForce();
	
	void ComputeForce();
	void ComputeForce(int i);
	//void DrawEdge(double a);

	void RK4Integration(double dt);
	void ImplicitIntegration(double dt);
	//void ExplicitEuler(double dt);

	void ConjugateGradientMethod(std::vector<std::vector<glm::mat3>> A, std::vector<glm::vec3> &v, std::vector<glm::vec3> b, int maxIterations, double tolerance);

	void BurstCoalesce(double dt);
	void UpdateAdjacentBubbles();
public:
	//solid
	void addSoild(glm::vec3 normal, glm::vec3 position);	
	void changeSolid(int index, glm::vec3 normal, glm::vec3 position);	
	void setParameter(float stiff, float weak, float dampv, float dampl, float solidad, float solidat, float burst, double wetness, float liquidAdhesion, float drag);

	//bubble type
	void setBubbleType(int i, float signed_distance);
	void setSurfaceNormal(int index, glm::vec3 surfaceNormal);
	void setLiquidVel(int i, glm::vec3 velocity);
	int getBubbleType(int i);
	void deleteBubble(int index);
	
//private:
	std::vector<CGAL::Triple<Weighted_point, Weighted_point, float>> alphaEdges;
	std::vector<CGAL::Triple<Weighted_point, Weighted_point, float>> non_alphaEdges;
	std::vector<Bubble> bubbleArr;
	std::vector<glm::vec3> solidNormals;
	std::vector<glm::vec3> solidPositions;
	std::vector<std::vector<int>> adjacency_matrix;
private:
	 float stiffCoeff;
	 float weakCoeff;
	 float dampCoeffV;
	 float dampCoeffL;
	 float solidAdhesionCoeff;
	 float solidAttractionCoeff;
	 float burstingCoeff;
	 float liquidAdhesionCoeff;
	 float dragCoeff;
	
	 double gl_wetness;
};
#endif



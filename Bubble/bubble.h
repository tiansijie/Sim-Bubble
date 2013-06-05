#ifndef BUBBLE_H_INCLUDED
#define BUBBLE_H_INCLUDED

//#include "openGL_headers.h"
#include "math_headers.h"
#include "glm.hpp"
#include <vector>
//#include "bubble.h"
#define SURFACE 0
#define LIQUID 1
#define AIR 2

class Bubble{
public:
	Bubble(){/*init_visualization();	*/gravity[1] = -9.8f;}
	Bubble(glm::vec3 center, glm::vec3 velocity, float radius):center(center), velocity(velocity), radius(radius){
		/*init_visualization();*/
		gravity[1] = -9.8f;
		weight = radius * radius;
		age = 0.0f;
		bubbleType = AIR;
		surfaceNormal = glm::vec3(0,1,0);
		mass = 4.0f/3.0f * 3.14f * pow(radius,3) * 1.29f;
		//mass = 1.0f;
		//mass = 1.0f;
		terminalSpeed = (2.0f / 3) * glm::sqrt(9.8 * radius);
		isSolidBubble = false;
	}
	glm::vec3 center;
	glm::vec3 velocity;
	float radius;
	float weight;
	float terminalSpeed;
	int degree;//use for Laplacian matrix
	float age;
	float burst;
	glm::vec3 ext_force;
	glm::vec3 gravity;
	glm::vec3 strong_interaction_force;
	glm::vec3 weak_interaction_force;
	glm::vec3 solid_attraction_force;
	glm::vec3 solid_adhesion_force;
	glm::vec3 damping_force;
	glm::vec3 buoyant_force;
	glm::vec3 drag_force;
	glm::vec3 liquid_adhesion_force;

	float signedDistance;
	int bubbleType;
	glm::vec3 surfaceNormal;
	glm::vec3 liquidVelocity;
	bool isSolidBubble;

	std::vector<int> adjacentBubbles;

	/*void draw(const VBO& vbo);*/
	//void update(glm::vec3 dis);
//private:
	/*void init_visualization();*/
	std::vector<glm::vec3> m_positions, m_colors, m_normals;
	std::vector<unsigned short> m_indices;
	float mass;
	
};
#endif
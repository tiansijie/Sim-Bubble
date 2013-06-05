#include "bubble_sim.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
//double gl_wetness = 0.0f;

static float mass = 1.0f;


std::vector<std::vector<glm::mat3>> M;//Mass matrix
glm::mat3 I = glm::mat3(1);//Identity matrix
std::vector<std::vector<glm::mat3>> J; //Jacobin matrix
std::vector<std::vector<glm::mat3>> J_DAMP;//Damping Jacobin matrix

Bubbles_sim::Bubbles_sim(){
	stiffCoeff = 2.0f;
	weakCoeff = 0.6f;
	dampCoeffV = 0.9f;
	dampCoeffL = 0.9f;
	solidAdhesionCoeff = 10.0f;
	solidAttractionCoeff = 6.0f;
	burstingCoeff = 0.02f;
	liquidAdhesionCoeff = 15.0f;
	dragCoeff = 0.5f;
	gl_wetness = 0.0f;	
}

Bubbles_sim::~Bubbles_sim(){
	
}

void Bubbles_sim::regularTriangulation_3()
{
	std::vector<Weighted_point> P;

//	int number_of_points = 0;
	int i = 0;
	for(std::vector<Bubble>::iterator it = bubbleArr.begin() ; it != bubbleArr.end(); ++it)
	{
		Point p(it->center.x, it->center.y, it->center.z);
		Weight w = it->weight ; 
		P.push_back(Weighted_point(p, w, i));
		i++;
//		++number_of_points;
	}

	Rt T;

	// insert all points in a row (this is faster than one insert() at a time).
	T.insert (P.begin(), P.end());

	assert( T.is_valid() );
	assert( T.dimension() == 3 );

	for(int i = 0; i < bubbleArr.size(); i++)
		bubbleArr[i].degree = 0;

	//initialize adjacency matrix
	adjacency_matrix.clear();
	adjacency_matrix.resize(bubbleArr.size());
	for(int i = 0; i < bubbleArr.size(); i++){
		bubbleArr[i].degree = 0;
		adjacency_matrix[i].clear();
		adjacency_matrix[i].resize(bubbleArr.size());
	}

//	std::cout << "Number of vertices : " << T.number_of_vertices() << std::endl;

	// removal of all vertices
//	int count = 0;

//	assert( count == number_of_points );
	Rt::Edge ed;
	Weighted_point p1, p2;
	
	alphaEdges.clear();
	non_alphaEdges.clear();

	for (Rt::Edge_iterator it = T.edges_begin(); it != T.edges_end(); ++it)
	{
		ed = *it;
		p1 = ed.first->vertex(ed.second)->point();
		p2 = ed.first->vertex(ed.third)->point();
		//std::cout<<p1.index()<<std::endl;
		int i1 = p1.index();
		int i2 = p2.index();
		if(i1 < 0 || i1 > bubbleArr.size() -1 || i2 < 0 || i2 > bubbleArr.size()-1)
			continue;
		bubbleArr[i1].degree++;//Degree for vertex
		bubbleArr[i2].degree++;//Degree for vertex

		//construct Adjacency matrix
		adjacency_matrix[i1][i2] = 1;
		adjacency_matrix[i2][i1] = 1;

		glm::vec3 vertex1(p1.x(), p1.y(), p1.z());
		glm::vec3 vertex2(p2.x(), p2.y(), p2.z());
		float radius1 = std::sqrt(p1.weight());
		float radius2 = std::sqrt(p2.weight());
		float restLen = std::sqrt(p1.weight() + p2.weight() + (3 * gl_wetness - 1) * radius1 * radius2);
		CGAL::Triple<Weighted_point, Weighted_point, float> edge(p1, p2, restLen);
		float distance = glm::distance(vertex1, vertex2);
		if(distance < (radius1 + radius2))
			alphaEdges.push_back(edge);
		else
			non_alphaEdges.push_back(edge);
	}

	
}

void Bubbles_sim::insertBubble(glm::vec3 position, glm::vec3 velocity, float radius){
	std::cout<<"b Position"<<position[0] << " "<< position[1] << " "<<position[2]<<std::endl;
	Bubble bu(position, velocity, radius);
	bubbleArr.push_back(bu);
}

Bubble Bubbles_sim::getBubble(int i){
	return bubbleArr.at(i);
}

glm::vec3 Bubbles_sim::getPosition(int i){
	return bubbleArr.at(i).center;
}

glm::vec3 Bubbles_sim::getVelocity(int i){
	return bubbleArr.at(i).velocity;
}

float Bubbles_sim::getRadius(int i){
	return bubbleArr.at(i).radius;
}

int Bubbles_sim::size(){
	return bubbleArr.size();
}

/*
void Bubbles_sim::draw(const VBO& vbos){
	for (int i = 0; i < bubbleArr.size(); i++)
	{
		bubbleArr[i].draw(vbos);
	}
}*/

/*
void Bubbles_sim::DrawEdge(double a){
	glBegin(GL_LINES);
	for(std::vector<CGAL::Triple<Weighted_point, Weighted_point, float>>::iterator it = non_alphaEdges.begin();
		it != non_alphaEdges.end(); it++){

			Weighted_point p1 = it->first;
			Weighted_point p2 = it->second;
			float rest_length = it->third;

			if(p1.index() < 0 || p2.index() < 0 || p1.index() >=  bubbleArr.size() || p2.index() >= bubbleArr.size())
				continue;

			glm::vec3 p1Vec3(p1.x(), p1.y(), p1.z());
			glm::vec3 p2Vec3(p2.x(), p2.y(), p2.z());

			glColor4f(1.0, 1.0, 0.0, a);
			//vec3 p1 = GetParticle(g, m_vsprings[i].m_p1).position;
			//vec3 p2 = GetParticle(g, m_vsprings[i].m_p2).position;
			glVertex3f(p1Vec3[0], p1Vec3[1], p1Vec3[2]);
			glVertex3f(p2Vec3[0], p2Vec3[1], p2Vec3[2]);
	}
	glEnd();
}*/

void Bubbles_sim::ComputeForce(){
//	DrawEdge(0.8f);
	SolidAdhesionForce();
	SolidAttractionForce();
	WeakBubbleForce();	
	BuoyantForce();
	DragForce();
	LiquidAdhesionForce();	
}

int coutttt = 0;

void Bubbles_sim::ComputeForce(int i){
}


void Bubbles_sim::update(double dt){	
	regularTriangulation_3();
	UpdateAdjacentBubbles();
	ImplicitIntegration(dt);
	//RK4Integration(dt);	
	BurstCoalesce(dt);
}


void Bubbles_sim::ImplicitIntegration(double dt){

	//computeForce();
	StrongBubbleForce();
	DampingForce();
	ComputeForce();

	//if(alphaEdges.size() != 0){

	float dt2 = dt*dt;
	std::vector<std::vector<glm::mat3>> A;
	//initialize A
	A.resize(bubbleArr.size());
	for(int i = 0; i < A.size(); i++)
		A[i].resize(bubbleArr.size());

	//Calculate A
	for(int i = 0; i < J.size(); i++){
		for(int j = 0; j < J.size(); j++){	
			//std::cout<< J[i][j][0][0] <<std::endl;
			A[j][i] = M[j][i] - (J[j][i] + J_DAMP[j][i]) * dt2;
		}
	}

	//Calculate b
	std::vector<glm::vec3> b;
	b.resize(bubbleArr.size());
	for(int i = 0; i < J.size(); i++){		
			b[i] = glm::vec3(0,0,0);		
	}

	for(int i = 0; i < J.size(); i++){
		//for(int j = 0; j < J.size(); j++){
		//	b[i] += (M[j][i] * bubbleArr[j].velocity);// + bubbleArr[i].strong_interaction_force * float(dt));////problem
		//}
		b[i] = M[i][i] * bubbleArr[i].velocity;
		/*std::cout<< "vel "<< bubbleArr[i].velocity[1] <<std::endl;
		std::cout<< "by "<< b[i][1] <<std::endl;*/
		//b[i][1] = - b[i][1];
		
		b[i] += (bubbleArr[i].strong_interaction_force + bubbleArr[i].damping_force
			+ bubbleArr[i].solid_attraction_force
			+ bubbleArr[i].solid_adhesion_force			
			+ bubbleArr[i].weak_interaction_force
			+ bubbleArr[i].buoyant_force
			+ bubbleArr[i].liquid_adhesion_force
			+ bubbleArr[i].drag_force
			+ bubbleArr[i].gravity * 4.0f/3.0f * 3.14f * pow(bubbleArr[i].radius,3) * 1.29f) * float(dt);	

		//std::cout<< "by "<< b[i][1] <<std::endl;
	}	
	std::vector<glm::vec3> vNew;
	vNew.resize(bubbleArr.size());

	ConjugateGradientMethod(A, vNew, b, 500, 0.0001);

	for(int i = 0; i < bubbleArr.size(); i++){	
		bubbleArr[i].center += (vNew[i] - bubbleArr[i].velocity) * float(dt);				
		bubbleArr[i].velocity = vNew[i];
	}
}

int notCoverTime = 0;

void Bubbles_sim::ConjugateGradientMethod(std::vector<std::vector<glm::mat3>> A, std::vector<glm::vec3> &v, std::vector<glm::vec3> b, int maxIterations, double tolerance){
	
	std::vector<glm::vec3> r;
	r.resize(bubbleArr.size());
	for(int i = 0; i < r.size(); i++){
		for(int j = 0; j < r.size(); j++){
			r[i] += A[j][i] * v[j];
		}		
		r[i] = b[i] - r[i];
	}
	
	std::vector<glm::vec3> d = r;
	std::vector<glm::vec3> q;
	
	float alpha_new = 0;
	float alpha = 0;
	float beta = 0;
	float delta_old = 0;
	float delta_new = 0;// = glm::dot(r, r);
	for(int i = 0; i < r.size(); i++){		
		delta_new += glm::dot(r[i], r[i]);
	}
	float delta0 = delta_new;

	float iter = 0;
	while(iter < maxIterations && delta_new > tolerance) {
		//std::cout<<"detlaNew is "<<delta_new<<std::endl;
		q.clear();
		q.resize(bubbleArr.size());
		//q = A*d;
		for(int i = 0; i < q.size(); i ++){
			for(int j = 0; j < q.size(); j++){
				q[i] += A[j][i] * d[j];
			}
		}

		//alpha = delta_new/dot(d,q);
		float temp = 0;
		for(int i = 0; i < d.size(); i++)
			temp += glm::dot(d[i],q[i]);
		alpha = delta_new/temp;

		//v = v + alpha*d;
		for(int i = 0; i < d.size(); i++)
			v[i] = v[i] + alpha * d[i];

		//r = r - alpha*q;
		for(int i = 0; i < r.size(); i++)
			r[i] = r[i] - alpha * q[i];

		delta_old = delta_new;

		//delta_new = glm::dot(r,r);
		delta_new = 0;
		for(int i = 0; i < r.size(); i++)
			delta_new += glm::dot(r[i], r[i]);

		beta = delta_new/delta_old;

		//d = r + beta*d;
		for(int i = 0; i < d.size(); i++)
			d[i] = r[i] + beta * d[i];
		iter++;

		/*if(delta_new < tolerance){
			std::cout<<"Less than tolerance"<<std::endl;
			break;
		}*/
	}	
	if(iter == maxIterations)
		std::cout<<"Wo diao"<<std::endl;
}

void Bubbles_sim::RK4Integration(double dt){
	//RK4 Method
	std::vector<glm::vec3> k1,k2,k3,k4;
	std::vector<glm::vec3> l1,l2,l3,l4;
	k1.resize(bubbleArr.size());
	k2.resize(bubbleArr.size());
	k3.resize(bubbleArr.size());
	k4.resize(bubbleArr.size());
	l1.resize(bubbleArr.size());
	l2.resize(bubbleArr.size());
	l3.resize(bubbleArr.size());
	l4.resize(bubbleArr.size());
	std::vector<glm::vec3> sV, sC;
	sV.reserve(bubbleArr.size());
	sC.reserve(bubbleArr.size());

	for(int i = 0; i< bubbleArr.size(); i++){
		sV[i] = bubbleArr[i].velocity;
		sC[i] = bubbleArr[i].center;
	}

	ComputeForce();
	for (int i = 0; i < bubbleArr.size(); i++)
	{
		ComputeForce(i);
		k1[i] = bubbleArr[i].velocity * float(dt);
		//std::cout<<"ext force "<<bubbleArr[i].ext_force[0] <<" "<< bubbleArr[i].ext_force[1] << " "<<bubbleArr[i].ext_force[2]<<std::endl;
		//std::cout<<" mass "<< bubbleArr[i].mass <<std::endl;
		l1[i] = bubbleArr[i].ext_force / bubbleArr[i].mass * float(dt);
		bubbleArr[i].velocity = sV[i] + l1[i] * 0.5f; 
		bubbleArr[i].center = sC[i] +  k1[i] * 0.5f;		
	}

	ComputeForce();
	for(int i = 0; i < bubbleArr.size(); i++){
		ComputeForce(i);
		k2[i] = bubbleArr[i].velocity * float(dt);// + l1[i] * float(dt) * 0.5f;
		l2[i] = bubbleArr[i].ext_force / bubbleArr[i].mass * float(dt);
		bubbleArr[i].velocity = sV[i] + l2[i] * 0.5f;
		bubbleArr[i].center = sC[i] +  k2[i] * 0.5f;		
	}

	ComputeForce();
	for(int i = 0; i < bubbleArr.size(); i++){
		ComputeForce(i);
		k3[i] = bubbleArr[i].velocity * float(dt);// + l2[i] * float(dt) * 0.5f;
		l3[i] = bubbleArr[i].ext_force / bubbleArr[i].mass * float(dt);
		bubbleArr[i].velocity = sV[i] + l3[i];
		bubbleArr[i].center = sC[i] +  k3[i];		
	}

	ComputeForce();
	for(int i = 0; i < bubbleArr.size(); i++){
		ComputeForce(i);
		k4[i] = bubbleArr[i].velocity * float(dt);// + l3[i] * float(dt);
		l4[i] = bubbleArr[i].ext_force / bubbleArr[i].mass * float(dt);
		//bubbleArr[i].velocity = sV[i] + l4[i] * 0.5f;
		//bubbleArr[i].center = sC[i] + k4[i];
	}

	for(int i = 0; i < bubbleArr.size(); i++){
		bubbleArr[i].velocity  += float(dt) / 6.0f * ( l1[i] + 2.0f * l2[i] + 2.0f * l3[i] + l4[i]);
		bubbleArr[i].center += float(dt) / 6.0f * ( k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]);
//		bubbleArr[i].update(float(dt) / 6.0f * ( k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]) + bubbleArr[i].center - sC[i]);
	}
}

//extern forces
void Bubbles_sim::StrongBubbleForce(){
	
	for(int i = 0; i < bubbleArr.size(); i++){
		bubbleArr[i].strong_interaction_force = glm::vec3(0,0,0);
	}

	//initialize for Jacobin and Mass matrix
	M.clear();
	M.resize(bubbleArr.size());
	J.clear();
	J.resize(bubbleArr.size());	
	for(int i = 0; i < J.size(); i++){
		M[i].resize(bubbleArr.size());
		J[i].resize(bubbleArr.size());
	}
	for(int i = 0; i< J.size(); i++)
		for(int j = 0; j< J.size(); j++)
			J[i][j] = glm::mat3(0);



	for(int i = 0; i < M.size(); i++){
		for(int j = 0 ; j < M.size(); j++){
			if(i != j)
				M[i][j] = glm::mat3(0);
			else if(i == j)
				M[i][j] = glm::mat3(bubbleArr[i].mass);
		}
	}

	
	//int i = 0;
	for(std::vector<CGAL::Triple<Weighted_point, Weighted_point, float>>::iterator it = alphaEdges.begin();
		it != alphaEdges.end(); it++){
			Weighted_point p1 = it->first;
			Weighted_point p2 = it->second;
			float rest_length = it->third;

			if(p1.index() < 0 || p2.index() < 0 || p1.index() >=  bubbleArr.size() || p2.index() >= bubbleArr.size())
				continue;

			glm::vec3 p1Vec3(p1.x(), p1.y(), p1.z());
			glm::vec3 p2Vec3(p2.x(), p2.y(), p2.z());

			//if(glm::length(p1Vec3-p2Vec3) <  rest_length)
				//std::cout<<"less than rest length"<<std::endl;

			bubbleArr[p1.index()].strong_interaction_force
				+=  -stiffCoeff * (p1Vec3 - p2Vec3 - rest_length * glm::normalize(p1Vec3 - p2Vec3));

			//bubbleArr[p1.index()].damping_force
				//+= (bubbleArr[p1.index()].velocity - bubbleArr[p2.index()].velocity);

			bubbleArr[p2.index()].strong_interaction_force
				+= -stiffCoeff * (p2Vec3 - p1Vec3 - rest_length * glm::normalize(p2Vec3 - p1Vec3));

			//bubbleArr[p2.index()].damping_force
			//	+= (bubbleArr[p2.index()].velocity - bubbleArr[p1.index()].velocity);

			//for the Jacobin matrix
			glm::vec3 deltP = p1Vec3 - p2Vec3;
			float dist = glm::length(deltP);
			if(dist >= rest_length){				
				int p1N = p1.index();
				int p2N = p2.index();
				J[p1N][p1N] -= (stiffCoeff * glm::outerProduct(deltP, deltP) / (dist*dist) + stiffCoeff * (1 - rest_length / dist) * (I - glm::outerProduct(deltP, deltP) / (dist * dist)));
				J[p2N][p1N] += (stiffCoeff * glm::outerProduct(deltP, deltP) / (dist*dist) + stiffCoeff * (1 - rest_length / dist) * (I - glm::outerProduct(deltP, deltP) / (dist * dist)));
				J[p2N][p2N] -= (stiffCoeff * glm::outerProduct(deltP, deltP) / (dist*dist) + stiffCoeff * (1 - rest_length / dist) * (I - glm::outerProduct(deltP, deltP) / (dist * dist)));
				J[p1N][p2N] += (stiffCoeff * glm::outerProduct(deltP, deltP) / (dist*dist) + stiffCoeff * (1 - rest_length / dist) * (I - glm::outerProduct(deltP, deltP) / (dist * dist)));
			}
			else{
				float k = (2 / rest_length) * ((dist/rest_length)/ sin(dist / rest_length));
				float fb = k * k * pow((cos(k*rest_length / 2.0f) - sin(k*rest_length/2.0f) / (k*rest_length/2.0f)),-1) * 1.0f/dist ;
				int p1N = p1.index();
				int p2N = p2.index();

				J[p1N][p1N] += 0;
				J[p2N][p1N] -= 0;
				J[p2N][p2N] += 0;
				J[p1N][p2N] -= 0;

				//deltP = p1Vec3 - p2Vec3;
				//dist = glm::length(deltP);
				
				/*if(fb < (stiffCoeff *(dist - rest_length))){
					J[p1N][p1N] += stiffCoeff * glm::outerProduct(deltP, deltP) / (dist*dist) +  (I - glm::outerProduct(deltP, deltP) / (dist * dist));
					J[p2N][p1N] += -(stiffCoeff * glm::outerProduct(deltP, deltP) / (dist*dist) + (I - glm::outerProduct(deltP, deltP) / (dist * dist)));
					J[p2N][p2N] += stiffCoeff * glm::outerProduct(deltP, deltP) / (dist*dist) + (I - glm::outerProduct(deltP, deltP) / (dist * dist));
					J[p1N][p2N] += -(stiffCoeff * glm::outerProduct(deltP, deltP) / (dist*dist) + (I - glm::outerProduct(deltP, deltP) / (dist * dist)));					
				}
				else{
					J[p1N][p1N] += fb * 1/dist * glm::outerProduct(deltP, deltP) / (dist*dist) +  (I - glm::outerProduct(deltP, deltP) / (dist * dist));
					J[p2N][p1N] += -(fb * 1/dist  * glm::outerProduct(deltP, deltP) / (dist*dist) + (I - glm::outerProduct(deltP, deltP) / (dist * dist)));
					J[p2N][p2N] += fb * 1/dist  * glm::outerProduct(deltP, deltP) / (dist*dist) + (I - glm::outerProduct(deltP, deltP) / (dist * dist));
					J[p1N][p2N] += -(fb * 1/dist  * glm::outerProduct(deltP, deltP) / (dist*dist) + (I - glm::outerProduct(deltP, deltP) / (dist * dist)));
				}*/
				//float fbStar = fb < (stiffCoeff *(dist - rest_length))? (stiffCoeff * (dist - rest_length)) : fb; 				
			}
	}
}

void Bubbles_sim::WeakBubbleForce(){
	
	for(int i = 0; i < bubbleArr.size(); i++)
		bubbleArr[i].weak_interaction_force = glm::vec3(0,0,0);

	//float coeff = 0.6f;

	for(std::vector<CGAL::Triple<Weighted_point, Weighted_point, float>>::iterator it = non_alphaEdges.begin();
		it != non_alphaEdges.end(); it++){

			Weighted_point p1 = it->first;
			Weighted_point p2 = it->second;
			float rest_length = it->third;

			if(p1.index() < 0 || p2.index() < 0 || p1.index() >=  bubbleArr.size() || p2.index() >= bubbleArr.size())
				continue;

			

			glm::vec3 p1Vec3(p1.x(), p1.y(), p1.z());
			glm::vec3 p2Vec3(p2.x(), p2.y(), p2.z());


			float sinij;
			float sinji;
			if(glm::distance(p1Vec3, p2Vec3) != 0.0f)
				sinij = bubbleArr[p2.index()].radius / glm::distance(p1Vec3, p2Vec3);
			else
				continue;
			if(glm::distance(p2Vec3, p1Vec3) != 0.0f)
				sinji = bubbleArr[p1.index()].radius / glm::distance(p2Vec3, p1Vec3);
			else
				continue;
			if(bubbleArr[p1.index()].bubbleType == SURFACE)
				bubbleArr[p1.index()].weak_interaction_force += weakCoeff * sinij * glm::normalize(p2Vec3 - p1Vec3);
			if(bubbleArr[p2.index()].bubbleType == SURFACE)
				bubbleArr[p2.index()].weak_interaction_force += weakCoeff * sinji * glm::normalize(p1Vec3 - p2Vec3);

			if(bubbleArr[p1.index()].weak_interaction_force[1] != bubbleArr[p1.index()].weak_interaction_force[1])
				bubbleArr[p1.index()].weak_interaction_force[1] = 0.0f;
			if(bubbleArr[p2.index()].weak_interaction_force[1] != bubbleArr[p2.index()].weak_interaction_force[1])
				bubbleArr[p2.index()].weak_interaction_force[1] = 0.0f;

	}
}

void Bubbles_sim::DragForce(){
	for(int i = 0; i < bubbleArr.size(); i++){
		if(bubbleArr[i].bubbleType == LIQUID)
		{
			bubbleArr[i].drag_force = dragCoeff * bubbleArr[i].radius * bubbleArr[i].radius * (bubbleArr[i].liquidVelocity - bubbleArr[i].velocity) * glm::length(bubbleArr[i].liquidVelocity - bubbleArr[i].velocity) / 100.0f ;
			//std::cout<<"Drag Force "<<bubbleArr[i].drag_force[0] << " "<< bubbleArr[i].drag_force[1]<<" " << bubbleArr[i].drag_force[2]<<std::endl;
		}
		else
			bubbleArr[i].drag_force = glm::vec3(0.0f, 0.0f, 0.0f);		
	}
}

void Bubbles_sim::BuoyantForce(){
	float liquidDensity = 1000.0f;
	for(int i = 0; i < bubbleArr.size(); i++){
		double cosV = (glm::dot(glm::normalize(bubbleArr[i].gravity), glm::normalize(bubbleArr[i].surfaceNormal)));		
		if(bubbleArr[i].bubbleType == LIQUID || (bubbleArr[i].bubbleType == SURFACE &&  cosV < 1.0f && cosV > -0.3f)){			
			bubbleArr[i].buoyant_force = - 4.0f/3.0f * 3.1415f * pow(bubbleArr[i].radius,3) * liquidDensity * bubbleArr[i].gravity;
		}
		else
			bubbleArr[i].buoyant_force = glm::vec3(0.0f, 0.0f, 0.0f);
	}
}

void Bubbles_sim::LiquidAdhesionForce(){
	for(int i = 0; i < bubbleArr.size(); i++){
		if(bubbleArr[i].bubbleType == SURFACE){ 
			bubbleArr[i].liquid_adhesion_force = -liquidAdhesionCoeff * bubbleArr[i].signedDistance * bubbleArr[i].surfaceNormal * bubbleArr[i].radius * 100.0f;//glm::vec3(0,1,0);
			//std::cout<<"Bubble Liquid " <<bubbleArr[i].liquid_adhesion_force[0] << " "<< bubbleArr[i].liquid_adhesion_force[1] << " "<< bubbleArr[i].liquid_adhesion_force[2]<<std::endl;
		}
		else
			bubbleArr[i].liquid_adhesion_force = glm::vec3(0,0,0);
	}
}

void Bubbles_sim::SolidAdhesionForce(){
	for(int i = 0; i < bubbleArr.size(); i++){
		bubbleArr[i].solid_adhesion_force = glm::vec3(0.0f, 0.0f, 0.0f);
		for(int j = 0; j < solidNormals.size(); j++){
			float signedDist = glm::dot(solidNormals[j], bubbleArr[i].center - solidPositions[j]);
			if(abs(signedDist) < bubbleArr[i].radius*2){ 
				bubbleArr[i].isSolidBubble = true;
				bubbleArr[i].solid_adhesion_force -= solidAdhesionCoeff * (signedDist - bubbleArr[i].radius*0.5f) * solidNormals[j];
			}
			else
				bubbleArr[i].isSolidBubble = false;
		}
	}
}

void Bubbles_sim::SolidAttractionForce(){
	for(int i = 0; i < bubbleArr.size(); i++){
		if(!bubbleArr[i].isSolidBubble){
			bubbleArr[i].solid_attraction_force = glm::vec3(0.0f, 0.0f, 0.0f);
			for(int j = 0; j < solidNormals.size(); j++){
				float dist = (glm::dot(solidNormals[j] , (bubbleArr[i].center - solidPositions[j])));
				//std::cout<<"distttt "<<dist<<std::endl;
				if(abs(dist) >  bubbleArr[i].radius)
					bubbleArr[i].solid_attraction_force -=  solidAttractionCoeff * solidNormals[j] / (dist - bubbleArr[i].radius*0.5f);
			}
		}
	}
}

void Bubbles_sim::DampingForce(){
	//construct Laplacian matrix
	std::vector<std::vector<float>> l_matrix;
	l_matrix.resize(bubbleArr.size());
	for(int i = 0; i < bubbleArr.size(); i++)
		l_matrix[i].resize(bubbleArr.size());

	for(int i = 0; i < bubbleArr.size(); i++){
		for(int j = 0; j < bubbleArr.size(); j++){
			if(i == j && bubbleArr[i].degree != 0)
				l_matrix[i][j] = 1.0f;
			else if(i != j && adjacency_matrix[i][j] == 1)
				l_matrix[i][j] = -1.0f / sqrt((float)bubbleArr[i].degree * bubbleArr[j].degree);
			else
				l_matrix[i][j] = 0.0f;
		}
	}
	
	std::vector<glm::vec3> term;
	term.resize(bubbleArr.size());
	for(int i = 0; i < bubbleArr.size(); i++){
		for(int j = 0; j < bubbleArr.size(); j++){
			term[i] += dampCoeffL * l_matrix[i][j] * bubbleArr[j].velocity;
		}
	}

	for(int i = 0; i < bubbleArr.size(); i++){
		bubbleArr[i].damping_force = -dampCoeffV * bubbleArr[i].velocity - term[i];
	}

	//initialize damping jacboin 
	J_DAMP.clear();
	J_DAMP.resize(bubbleArr.size());
	for(int i = 0; i < bubbleArr.size(); i++){
		J_DAMP[i].clear();
		J_DAMP[i].resize(bubbleArr.size());
	}

	for(int i = 0; i < bubbleArr.size(); i++)
		for(int j = 0; j < bubbleArr.size(); j++)
			J_DAMP[j][i] = glm::mat3(0);

	for(int i = 0; i < bubbleArr.size(); i++){
		for(int j = 0; j < bubbleArr.size(); j++){
			if(i == j)
				J_DAMP[j][i] = -(dampCoeffV + dampCoeffL) * I;
			else if(i != j && adjacency_matrix[i][j] == 1)
				J_DAMP[j][i] = -dampCoeffL * 1 / sqrt((float)bubbleArr[i].degree * bubbleArr[j].degree) * I;
			else
				J_DAMP[j][i] = glm::mat3(0);
		}
	}
}


void Bubbles_sim::addSoild(glm::vec3 normal, glm::vec3 position){
	solidNormals.push_back(normal);
	solidPositions.push_back(position);
}


void Bubbles_sim::changeSolid(int index, glm::vec3 normal, glm::vec3 position){
	solidNormals[index] = normal;
	solidPositions[index] = position;
}

void Bubbles_sim::setParameter(float stiff, float weak, float dampv, float dampl, float solidad, float solidat, float burst, double wetness, float liquidAdhesion, float drag){
	stiffCoeff = stiff; weakCoeff = weak; dampCoeffV = dampv; dampCoeffL = dampl; solidAdhesionCoeff = solidad;
	solidAttractionCoeff = solidat; burstingCoeff = burst; gl_wetness = wetness; liquidAdhesionCoeff = liquidAdhesion; dragCoeff = drag;
}


void Bubbles_sim::BurstCoalesce(double dt){
	for(int i = 0; i < bubbleArr.size(); i++){
		bubbleArr[i].age += dt*2.0f;
		bubbleArr[i].burst = 1 - 1 / exp(burstingCoeff * bubbleArr[i].radius * bubbleArr[i].age);
		float minDistance = DBL_MAX;
		int minBubble = -1;
		float kkk = 1.0f - bubbleArr[i].burst;		
		if(1.0f - bubbleArr[i].burst < 0.01f + (rand() %100)/1000.0f){
		//	std::cout<<"It's bursting !!!!!!!!!"<<std::endl;
			for(int j = 0; j < bubbleArr[i].adjacentBubbles.size(); j++){//find min distance
				int n = bubbleArr[i].adjacentBubbles[j];
				float dist = glm::length(bubbleArr[i].center - bubbleArr[n].center);
				if(dist < minDistance){
					minDistance = dist;
					minBubble = n;
				}
			}
			if(minBubble != -1){
				bubbleArr[minBubble].center = (0.3f * bubbleArr[i].center + 0.7f *bubbleArr[minBubble].center);
				bubbleArr[minBubble].radius = pow(pow(bubbleArr[minBubble].radius,3) + pow(bubbleArr[i].radius/2.0f,3),1.0f/3.0f);				
			}
			bubbleArr.erase(bubbleArr.begin() + i);
		}
	}
}

void Bubbles_sim::UpdateAdjacentBubbles(){
	for(int i = 0; i < bubbleArr.size(); i++){
		bubbleArr[i].adjacentBubbles.clear();
	}

	for(std::vector<CGAL::Triple<Weighted_point, Weighted_point, float>>::iterator it = alphaEdges.begin();it != alphaEdges.end(); it++){

			Weighted_point p1 = it->first;
			Weighted_point p2 = it->second;
			float rest_length = it->third;

			if(p1.index() < 0 || p2.index() < 0 || p1.index() >=  bubbleArr.size() || p2.index() >= bubbleArr.size())
				continue;

			bubbleArr[p1.index()].adjacentBubbles.push_back(p2.index());
			bubbleArr[p2.index()].adjacentBubbles.push_back(p1.index());
	}
}

void Bubbles_sim::setBubbleType(int i, float signed_distance){
	bubbleArr[i].signedDistance = signed_distance;
	if(signed_distance < -bubbleArr[i].radius){
		std::cout<<"Liquid bubble "<<std::endl;
		bubbleArr[i].bubbleType = LIQUID;
	}
	else if(signed_distance > -bubbleArr[i].radius && signed_distance < bubbleArr[i].radius){
		std::cout<<"Surface bubble "<<std::endl;
		bubbleArr[i].bubbleType = SURFACE;
	}
	else{
		std::cout<<"Air bubble "<<std::endl;
		bubbleArr[i].bubbleType = AIR;
	}
}

void Bubbles_sim::setSurfaceNormal(int index, glm::vec3 normal){
	bubbleArr[index].surfaceNormal = glm::normalize(normal);
}

void Bubbles_sim::setLiquidVel(int i, glm::vec3 velocity){
	bubbleArr[i].liquidVelocity = velocity;
}


void Bubbles_sim::deleteBubble(int index){
	bubbleArr.erase(bubbleArr.begin() + index);
}

int Bubbles_sim::getBubbleType(int i){
	return bubbleArr[i].bubbleType;
}
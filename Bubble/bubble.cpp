#include "bubble.h"


/*
void Bubble::update(glm::vec3 dis){
	/ *for(int i = 0; i < bubbleArr.size(); i++){* /
	//if((center.y) > 0.5f){
	//	ext_force = gravity + solid_attraction_force + solid_adhesion_force + strong_interaction_force + weak_interaction_force;
	//}
	//else{
	//	ext_force = solid_attraction_force + solid_adhesion_force + strong_interaction_force + weak_interaction_force - damping_force - velocity; 
	//	velocity.y = 0.0f;
	//}
	//	//glm::vec3 tp, tv;

	//	//tv = velocity + ext_force * (float) dt * 0.5f;
	//	//tp = center + velocity * (float) dt * 0.5f;


	//	glm::vec3 k1,k2,k3,k4;
	//	glm::vec3 l1,l2,l3,l4;

	//	k1 = velocity;
	//	
	//	k2 = velocity + ext_force * float(dt) * 0.5f;

	//	k3 = velocity + ext_force * float(dt) * 0.5f;

	//	k4 = velocity + ext_force * float(dt);
	//	
	//	velocity  += ext_force * (float)dt;//+= ext_force * (float)dt;

	//	
	//	for(int i = 0; i < m_positions.size(); i++){
	//		m_positions[i] += float(dt) / 6.0f * ( k1 + 2.0f * k2 + 2.0f * k3 + k4);//tv * float(dt);//(velocity + k2) / 2.0f;//velocity  * (float)dt;
	//	}
	//	
	//	center += float(dt) / 6.0f * ( k1 + 2.0f * k2 + 2.0f * k3 + k4);//tv * float(dt);//(velocity + k2) / 2.0f;//velocity * (float)dt;
	////}	 
	
	for(int i = 0; i < m_positions.size(); i++){
			m_positions[i] += dis;//tv * float(dt);//(velocity + k2) / 2.0f;//velocity  * (float)dt;
	}


}*/

/*
void Bubble::init_visualization()
{
	m_positions.clear();
	m_colors.clear();
	m_normals.clear();
	m_indices.clear();

	glm::vec3 mat_color(0.2f, 0.6f, 0.5f);
	unsigned int slice = 24, stack = 10;


	glm::vec3 tnormal(0.0f, 1.0f, 0.0f), tpos;
	//tpos = m_center + m_radius * tnormal;
	tpos = center + radius * tnormal;

	m_positions.push_back(tpos);
	m_normals.push_back(tnormal);
	m_colors.push_back(mat_color);

	float theta_z, theta_y, sin_z;
	float delta_y = 360.0f / slice, delta_z = 180.0f / stack;
	//loop over the sphere
	for(theta_z = delta_z; theta_z < 179.99f; theta_z += delta_z)
	{
		for(theta_y = 0.0f; theta_y < 359.99f; theta_y += delta_y)
		{
			sin_z = sin(glm::radians(theta_z));

			tnormal.x = sin_z * cos(glm::radians(theta_y));
			tnormal.y = cos(glm::radians(theta_z));
			tnormal.z = -sin_z * sin(glm::radians(theta_y));

			tpos = center + radius * tnormal;

			m_positions.push_back(tpos);
			m_normals.push_back(tnormal);
			m_colors.push_back(mat_color);
		}
	}
	tnormal = glm::vec3(0.0f, -1.0f, 0.0f);
	tpos = center + radius * tnormal;

	m_positions.push_back(tpos);
	m_normals.push_back(tnormal);
	m_colors.push_back(mat_color);

	//indices
	unsigned int j = 0, k = 0;
	for(j = 0; j < slice - 1; ++j)
	{
		m_indices.push_back(0);
		m_indices.push_back(j + 1);
		m_indices.push_back(j + 2);
	}
	m_indices.push_back(0);
	m_indices.push_back(slice);
	m_indices.push_back(1);

	for(j = 0; j < stack - 2; ++j)
	{
		for(k = 1 + slice * j; k < slice * (j + 1); ++k)
		{
			m_indices.push_back(k);
			m_indices.push_back(k + slice);
			m_indices.push_back(k + slice + 1);

			m_indices.push_back(k);
			m_indices.push_back(k + slice + 1);
			m_indices.push_back(k + 1);
		}
		m_indices.push_back(k);
		m_indices.push_back(k + slice);
		m_indices.push_back(k + 1);

		m_indices.push_back(k);
		m_indices.push_back(k + 1);
		m_indices.push_back(k + 1 - slice);
	}

	unsigned int bottom_id = (stack - 1) * slice + 1;
	unsigned int offset = bottom_id - slice;
	for(j = 0; j < slice - 1; ++j)
	{
		m_indices.push_back(j + offset);
		m_indices.push_back(bottom_id);
		m_indices.push_back(j + offset + 1);
	}
	m_indices.push_back(bottom_id - 1);
	m_indices.push_back(bottom_id);
	m_indices.push_back(offset);

	if(m_indices.size() != 6 * (stack - 1) * slice)
		printf("indices number not correct!\n");

}

void Bubble::draw(const VBO& vbos)
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	// position
	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
	glBufferData(GL_ARRAY_BUFFER, 3 * m_positions.size() * sizeof(float), &m_positions[0], GL_STREAM_DRAW);

	// color
	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
	glBufferData(GL_ARRAY_BUFFER, 3 * m_colors.size() * sizeof(float), &m_colors[0], GL_STREAM_DRAW);

	// normal
	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
	glBufferData(GL_ARRAY_BUFFER, 3 * m_normals.size() * sizeof(float), &m_normals[0], GL_STREAM_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned short), &m_indices[0], GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);

	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
	glDrawElements(GL_TRIANGLES, m_indices.size(), GL_UNSIGNED_SHORT, 0);//GL_UNSIGNED_INT

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}*/
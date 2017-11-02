/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <map>
#include <algorithm> //Para fill_n
#include <vector>

#include "particle_filter.h"
#include <math.h>


using namespace std;

Particle part;
default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).


	
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	num_particles = 200;

	for (int i = 0; i<num_particles; i++)
	{
		//double sample_x, sample_y, sample_theta;
		//string name = "part";
		//+ std::to_string(i);
		//Particle particles;
		Particle part;
		part.id = num_particles;
		part.x = dist_x(gen);
		part.y = dist_y(gen);
		part.theta = dist_theta(gen);
		particles.push_back(part);

		//cout << "Sample" << " " << part.x  << " " << part.y  << " " << part.theta  << endl;

	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	
	normal_distribution<double> Noise_x(0, std_pos[0]);
	normal_distribution<double> Noise_y(0, std_pos[1]);
	normal_distribution<double> Noise_theta(0, std_pos[2]);

	for (int i=0; i <num_particles; i++)
	{
		if (fabs(yaw_rate)<0.00001){
			particles[i].x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y = particles[i].y + velocity/yaw_rate*(-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta));
		}
		else {
			particles[i].x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y = particles[i].y + velocity/yaw_rate*(-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta));
			particles[i].theta = particles[i].theta + yaw_rate*delta_t;
		}

		//Add noise
		particles[i].x += Noise_x(gen);
		particles[i].y += Noise_y(gen);
		particles[i].theta += Noise_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


	//for each observation it compares with each predicted and calculates the minimum distance
	
	for (auto& observation:observations){
			double min_distance = 99999999.0;
			for (const auto& p:predicted) {
				double distance = dist(p.x, p.y, observation.x, observation.y);
				if(distance < min_distance){
					min_distance = distance;
					observation.id = p.id;
				}
			}
	} 
	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) { 
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	//Transformation of the car observations
	
	double x_map;
	double y_map;

	vector<LandmarkObs> transformed_obs; //vehicle coordinates
	transformed_obs.resize(observations.size()); //vehicle coordinates?

	for (int i=0; i < num_particles; i++){
	Particle &particle = particles[i];
	

		for (int j=0; j<observations.size(); j++){

			LandmarkObs transformed_ob;

			transformed_ob.x = particle.x + (cos(particle.theta) * observations[j].x) - (sin(particle.theta)* observations[j].y);
			transformed_ob.y = particle.y + (sin(particle.theta) * observations[j].x) + (cos(particle.theta)* observations[j].y);
			
			transformed_obs[j] = transformed_ob;
		}

		std::vector<LandmarkObs> inrange_landmarks; 
		std::map<int, Map::single_landmark_s> idx2landmark;

		for (const auto& landmark:map_landmarks.landmark_list){ //I grab a landmark from landmark_list
			double distance = dist(landmark.x_f, landmark.y_f, particle.x, particle.y);
			if(distance <= sensor_range){
				inrange_landmarks.push_back(LandmarkObs{landmark.id_i, landmark.x_f, landmark.y_f});
				idx2landmark.insert(std::make_pair(landmark.id_i, landmark)); //Set x and y
			}
		}

		if (inrange_landmarks.size() > 0) {

			dataAssociation(inrange_landmarks, transformed_obs);

			particle.weight = 1.0; //reset the weight of the particle

		



	//Calculate the weight of each particle
	float gauss_norm;
	float exponent;
	float weight;

	for (const auto observation:transformed_obs){
		gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
		exponent = (pow(observation.x - idx2landmark[observation.id].x_f, 2))/(2 * pow(std_landmark[0],2)) + (pow(observation.y - idx2landmark[observation.id].y_f, 2))/(2 * pow(std_landmark[1], 2));
		particle.weight *= gauss_norm * exp(-exponent); 
		//cout << particle.weight << endl;
	}
	//Segmentation fault weights[i] = particle.weight;
	weights.push_back(particle.weight); 
	//cout << weights[i] << endl; //Imprime el peso de cada particula

	
}
else{
	//Segmentation fault weights[i] = 0.0;
	weights.push_back(0.0);
	//cout << weights[i] << endl; //Imprime cuando el peso es cero
}
}
		} 


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	
	default_random_engine gen;
	discrete_distribution<int> dist_p(weights.begin(), weights.end()); //como enumerar las particulas
	//fill(weights.begin(),weights.end(),0);
	//std::fill_n(weights.bi, 200, 0); //Error con los operadores
	

	vector<Particle> new_particles;

	//Size of weight
	cout << "Tamaño de weights: "<< weights.size() << endl;

	//Size of particles
	cout << "Tamaño de particles: " << particles.size() << endl;

	weights.resize(0);

	for(int i=0; i<num_particles; i++){
		new_particles.push_back(particles[dist_p(gen)]);
		new_particles[i].id = i;
	}
	weights.resize(0);
	particles = new_particles;
	
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

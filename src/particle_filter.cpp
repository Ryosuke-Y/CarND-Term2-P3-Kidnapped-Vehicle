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

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 200;

  particles.resize(num_particles);
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0,1);

  for(int i=0; i<num_particles; i++){

    particles[i].x = x + std[0] * distribution(generator);
    particles[i].y = y + std[1] * distribution(generator);
    particles[i].theta = theta + std[2] * distribution(generator);
    particles[i].weight = 1.;
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std::default_random_engine gen;

  for(int i=0;i<num_particles;i++){

    Particle &p = particles[i];

    double x_pred;
    double y_pred;
    double theta_pred;

    theta_pred = p.theta + yaw_rate*delta_t;
    x_pred = p.x + velocity/yaw_rate*(sin(theta_pred)-sin(p.theta));
    y_pred = p.y - velocity/yaw_rate*(cos(theta_pred)-cos(p.theta));

    std::normal_distribution<double> dist_x(x_pred,std_pos[0]);
    std::normal_distribution<double> dist_y(y_pred,std_pos[1]);
    std::normal_distribution<double> dist_theta(theta_pred,std_pos[2]);

    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);

  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for (unsigned i = 0; i < observations.size(); ++i) {

    // Define temporary variables for finding predicted measurements
    double dist_current = 1e6;
    int nearest_landmark = -1;

    // Iterate through predicted measurements to find nearest landmarks
    for (unsigned j = 0; j < predicted.size(); ++j) {

      // Calculate Euclidian distance between observations and predictions
      double dist_eucl = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

      // assign min
      if (dist_eucl < dist_current) {
        dist_current = dist_eucl;
        nearest_landmark = j;
      }
    }
    // assign the closest id to the obeservation
    observations[i].id = predicted[nearest_landmark].id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::default_random_engine generator;
	std::discrete_distribution<int> dist(weights.begin(), weights.end());

	std::vector<Particle> new_particles;
	new_particles.reserve(particles.size());

	for (int i = 0; i < particles.size(); i++) {
		int sampled_index = dist(generator);
		new_particles.push_back(particles[sampled_index]);
	}

	particles = new_particles;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}

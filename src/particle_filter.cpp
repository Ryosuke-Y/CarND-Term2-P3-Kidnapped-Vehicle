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
	num_particles = 100;

  particles.resize(num_particles);
	weights.resize(num_particles);

  std::default_random_engine generator;
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

  for(int i=0; i<num_particles; i++){

		particles[i].id = i;
		particles[i].x = dist_x(generator);
		particles[i].y = dist_y(generator);
		particles[i].theta = dist_theta(generator);
		particles[i].weight = 1.0;
		weights[i] = 1.0;
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std::default_random_engine generator;

	double x_pred;
	double y_pred;
	double theta_pred;

  for(int i=0;i<num_particles;i++){

    Particle &p = particles[i];

		if (fabs(yaw_rate) > 0.0001) {

			x_pred = p.x + velocity/yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
			y_pred = p.y + velocity/yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
			theta_pred = p.theta + yaw_rate * delta_t;

		} else {

			x_pred = p.x + velocity * delta_t * cos(p.theta);
			y_pred = p.y + velocity * delta_t * sin(p.theta);
			theta_pred = p.theta;

		}

		// Update the particle position with the prediction and add gaussian noise
		std::normal_distribution<double> dist_x(x_pred, std_pos[0]);
		std::normal_distribution<double> dist_y(y_pred, std_pos[1]);
		std::normal_distribution<double> dist_theta(theta_pred, std_pos[2]);

		p.x = dist_x(generator);
		p.y = dist_y(generator);
		p.theta = dist_theta(generator);

  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i=0; i<observations.size(); ++i) {

        double dist_min = 0.0;
        int id_m = 0;

        for (int j=0; j<predicted.size(); ++j) {

            double dist_temp = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

            if (j == 0) {
                dist_min = dist_temp;
                id_m = predicted[j].id;
            }

            if (dist_temp < dist_min) {
                dist_min = dist_temp;
                id_m = predicted[j].id;
            }
        }

        observations[i].id = id_m;
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

	double weights_sum = 0.0;

	for (int i = 0; i < particles.size(); i++) {
		Particle &p = particles[i];

		// Make a list of predicted measurements for this particle
		std::vector<LandmarkObs> predicted(map_landmarks.landmark_list.size());
		for (int j = 0; j < predicted.size(); j++) {
			Map::single_landmark_s landmark = map_landmarks.landmark_list[j];

			double sin_theta = sin(p.theta);
			double cos_theta = cos(p.theta);

			double x_pred = (landmark.y_f- p.y)*sin_theta + (landmark.x_f - p.x) * cos_theta;
			double y_pred = (landmark.y_f- p.y)*cos_theta - (landmark.x_f - p.x) * sin_theta;

			predicted[j].x = x_pred;
			predicted[j].y = y_pred;
			predicted[j].id = landmark.id_i;
		}

		dataAssociation(predicted, observations);

		double prob = 1;

		for (LandmarkObs &obs : observations) {

			double obs_x = obs.x * cos(p.theta) - obs.y * sin(p.theta) + p.x;
			double obs_y = obs.x * sin(p.theta) + obs.y * cos(p.theta) + p.y;

			double true_x = map_landmarks.landmark_list[obs.id-1].x_f;
			double true_y = map_landmarks.landmark_list[obs.id-1].y_f;

			double obs_prob = (1/(2*M_PI*std_landmark[0]*std_landmark[1]))
				* exp(-(pow(true_x - obs_x, 2)/(2*std_landmark[0]*std_landmark[0])
					+ (pow(true_y - obs_y, 2)/(2*std_landmark[1]*std_landmark[1]))));

			prob *= obs_prob;
		}

		weights[i] = prob;
	}

  // normalize weights
  for (int i=0; i<num_particles; i++)
      particles[i].weight /= weights_sum;

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
		Particle &p = particles[i];
		dataFile << p.x << " " << p.y << " " << p.theta << "\n";
	}
	dataFile.close();
}

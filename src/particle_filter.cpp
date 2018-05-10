/*
 * particle_filter.cpp
 *
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#define THRESHOLD 0.0001

#include "particle_filter.h"

using namespace std;

void ParticleFilter::Init(const double & x, const double & y, const double & theta, const double *std) {
  // Set the number of particles. Initialize all particles to first position (based on estimates of
  // x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  ///* Tuning parameters

  num_particles = 101;

  ///* Tuning parameters

  // Resize weights vector
  weights.resize(num_particles);

  // Resize particles vector
  particles.resize(num_particles);

  // Engine for generating the particles randomly
  random_device random_device_;
  default_random_engine random_engine_(random_device_());

  // Create normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> norm_dist_x(x, std[0]);
  normal_distribution<double> norm_dist_y(y, std[1]);
  normal_distribution<double> norm_dist_theta(theta, std[2]);

  // Initialize the particles
  int _id = 0;
  for (auto &particle : particles) {
    particle.id = _id++;
    particle.x = norm_dist_x(random_engine_);
    particle.y = norm_dist_y(random_engine_);
    particle.theta = norm_dist_theta(random_engine_);
    particle.weight = 1.0;

    weights[_id - 1] = particle.weight;
  }

  // Done!
  is_initialized = true;
}

void ParticleFilter::Prediction(const double &delta_t, const double *std_pos, const double &velocity, const double &yaw_rate) {
  // Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  // * http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  // * http://www.cplusplus.com/reference/random/default_random_engine/

  // Engine for generating the particles randomly
  random_device random_device_;
  default_random_engine random_engine_(random_device_());

  // Create normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> norm_dist_x(0.0, std_pos[0]);
  normal_distribution<double> norm_dist_y(0.0, std_pos[1]);
  normal_distribution<double> norm_dist_theta(0.0, std_pos[2]);

  // To have less calculation and improve efficiency
  const auto vdt = velocity * delta_t;
  const auto yawddt = yaw_rate * delta_t;
  auto v_yawd = 0.0;
  if (fabs(yaw_rate) >= THRESHOLD) {
    v_yawd = velocity / yaw_rate;
  }

  // Add measurements to particles
  for (auto &particle : particles) {
    // Avoid dividing by ZERO
    if (fabs(yaw_rate) < THRESHOLD) {
      particle.x += vdt * cos(particle.theta);
      particle.y += vdt * sin(particle.theta);
    } else {
      const auto arg = particle.theta + yawddt;
      particle.x += v_yawd * (sin(arg) - sin(particle.theta));
      particle.y += v_yawd * (-cos(arg) + cos(particle.theta));
      particle.theta += yawddt;
    }

    // Add normal noise
    particle.x += norm_dist_x(random_engine_);
    particle.y += norm_dist_y(random_engine_);
    particle.theta += norm_dist_theta(random_engine_);
  }
}

void ParticleFilter::DataAssociation(const std::vector<LandmarkObs> &predicted,
                                     std::vector<LandmarkObs> &observations) {
  // Find the predicted measurement that is closest to each observed measurement and assign the
  // observed measurement to this particular landmark.
  // NOTE: This method will NOT be called by the grading code. But you will probably find it useful to
  // implement this method and use it as a helper during the UpdateWeights phase.

  for (auto &obs : observations) {
    // Set minimum distance to maximum possible value fro "double"
    double min_dist = numeric_limits<double>::max();

    // Init id of landmark from map placeholder to be associated with the observation
    int nearest_lm_id = -1;

    for (auto &pred : predicted) {
      // Get distance between current/predicted landmarks
      double cur_dist = dist(obs.x, obs.y, pred.x, pred.y);

      // Find the predicted landmark nearest the current observed landmark
      if (cur_dist < min_dist) {
        min_dist = cur_dist;
        nearest_lm_id = pred.id;
      }
    }

    // Associate the observation's id with the nearest predicted landmark's identification
    obs.id = nearest_lm_id;
  }

}

void ParticleFilter::TransformToMapCoordinate(const Particle &particle, const std::vector<LandmarkObs> &observations,
                                              std::vector<LandmarkObs> &transformed_obs) {
  for (auto &obs : observations) {
    // Homogeneous Transformation
    const double transformed_obs_x = cos(particle.theta) * obs.x - sin(particle.theta) * obs.y + particle.x;
    const double transformed_obs_y = sin(particle.theta) * obs.x + cos(particle.theta) * obs.y + particle.y;
    transformed_obs.push_back(LandmarkObs{obs.id, transformed_obs_x, transformed_obs_y});
  }
}

void ParticleFilter::FindInRangeLandmarks(const double &sensor_range, const Particle &particle,
                                          const Map &map_landmarks, std::vector<LandmarkObs> &in_range_particles) {
  for (auto &lm : map_landmarks.landmark_list) {
    const auto dist_ = dist(lm.x_f, lm.y_f, particle.x, particle.y);
    if (dist_ <= sensor_range) {
      in_range_particles.push_back(LandmarkObs{lm.id_i, lm.x_f, lm.y_f});
    }
  }
}

double ParticleFilter::CalculateParticleWeight(const double *std_landmark,
                                               const std::vector<LandmarkObs> &predicted,
                                               const std::vector<LandmarkObs> &observations) {

  const double std_x = std_landmark[0];
  const double std_y = std_landmark[1];

  // Initialize the weight
  double weight = 1.0;

  for (auto &obs : observations) {
    const int obs_id = obs.id;

    for (auto &pred : predicted) {
      if (pred.id == obs_id) {
        // Final weight is calculated by multiplying all the calculated measurement probabilities together
        weight *= this->MultivariateGaussianProbability(obs.x, obs.y, pred.x, pred.y, std_x, std_y);
        break;
      }
    }
  }

  return weight;
}

double ParticleFilter::MultivariateGaussianProbability(const double &x, const double &y, const double &mu_x,
                                                       const double &mu_y, const double &std_x, const double &std_y) {
  return 1.0 / (2.0 * M_PI * std_x * std_y) *
         exp(-0.5 * pow(x - mu_x, 2) / pow(std_x, 2) - 0.5 * pow(y - mu_y, 2) / pow(std_y, 2));
}

void ParticleFilter::UpdateWeights(const double &sensor_range, const double *std_landmark,
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
  // Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  // more about this distribution here:
  // * https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  //
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  // according to the MAP'S coordinate system. You will need to transform between the two systems.
  // Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //
  // The following is a good resource for the theory:
  // * https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  // and the following is a good resource for the actual equation to implement (look at equation 3.33):
  // * http://planning.cs.uiuc.edu/node99.html

  // This variable is used for normalizing weights of all particles and bring them in the range
  // of [0, 1]
  double weight_normalizer = 0.0;

  for (auto &particle : particles) {
    // Transform landmark objects from the vehicle coordinate to the map coordinate
    vector<LandmarkObs> transformed_obs;
    this->TransformToMapCoordinate(particle, observations, transformed_obs);

    // Define a holder for landmarks which are in sensor range
    vector<LandmarkObs> predictions;
    this->FindInRangeLandmarks(sensor_range, particle, map_landmarks, predictions);

    // Associate data for the predictions and transformed observations on the particle
    this->DataAssociation(predictions, transformed_obs);

    // Calculate the particle's weight
    particle.weight = this->CalculateParticleWeight(std_landmark, predictions, transformed_obs);

    weight_normalizer += particle.weight;
  }

  // Normalize the weights of all particles since re-sampling using probabilistic approach.
  unsigned int idx = 0;
  if (weight_normalizer > 0) {
    //cout << ">> weight_normalizer: " << weight_normalizer << endl;
    for (auto &particle: particles) {
      particle.weight /= weight_normalizer;
      weights[idx++] = particle.weight;
    }
  }
}

void ParticleFilter::Resample() {
  // Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  // * http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // Engine for generating the particles randomly
  random_device random_device_;
  default_random_engine random_engine_(random_device_());
  // Generate random starting index for resampling wheel
  uniform_int_distribution<long> uni_int_dist(0, num_particles - 1);
  auto index = uni_int_dist(random_engine_);

  // Get max weight
  double max_weight_2 = *max_element(weights.begin(), weights.end()) * 2.0;
  // Uniform random distribution [0.0, max_weight)
  uniform_real_distribution<double> uni_real_dist(0.0, max_weight_2);

  vector<Particle> re_sampled_particles;
  double beta = 0.0;

  // Spin the resample wheel!
  auto num_p = static_cast<unsigned int>(num_particles);
  for (unsigned int i = 0; i < num_p; i++) {
  //for (auto &_ : weights) {
    beta += uni_real_dist(random_engine_);
    while (beta > weights[index]) {
      beta -= weights[index];
      index = static_cast<int>((index + 1) % num_particles);
    }
    re_sampled_particles.push_back(particles[index]);
  }

  particles = re_sampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle &particle, const std::vector<int> &associations,
                                         const std::vector<double> &sense_x, const std::vector<double> &sense_y) {
  // particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;

  return particle;
}

string ParticleFilter::GetAssociations(Particle best) {
  vector<int> v = best.associations;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::GetSenseX(Particle best) {
  vector<double> v = best.sense_x;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::GetSenseY(Particle best) {
  vector<double> v = best.sense_y;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
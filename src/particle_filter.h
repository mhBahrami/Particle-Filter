/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

struct Particle {

  int id;
  double x;
  double y;
  double theta;
  double weight;

  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;

};


class ParticleFilter {

  // Number of particles to draw
  unsigned long num_particles;


  // Flag, if filter is Initialized
  bool is_initialized;

  // Vector of weights of all particles
  std::vector<double> weights;

public:

  /**
   * Set of current particles
   */
  std::vector<Particle> particles;

  /**
   * Constructor
   * @param num_particles Number of particles
   */
  ParticleFilter() : num_particles(0), is_initialized(false) {}

  /**
   * Destructor
   */
  ~ParticleFilter() = default;

  /**
   * Init Initializes particle filter by initializing particles to Gaussian
   * distribution around first position and all the weights to 1.
   * @param x Initial x position [m] (simulated estimate from GPS)
   * @param y Initial y position [m]
   * @param theta Initial orientation [rad]
   * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
   *   standard deviation of yaw [rad]]
   */
  void Init(const double &x, const double &y, const double &theta, const double *std);

  /**
   * Prediction Predicts the state for the next time step
   * using the process model.
   * @param delta_t Time between time step t and t+1 in measurements [s]
   * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
   *        standard deviation of yaw [rad]]
   * @param velocity Velocity of car from t to t+1 [m/s]
   * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
   */
  void Prediction(const double &delta_t, const double *std_pos, const double &velocity, const double &yaw_rate);

  /**
   * DataAssociation Finds which observations correspond to which landmarks (likely by using
   * a nearest-neighbors data association).
   * @param predicted Vector of predicted landmark observations
   * @param observations Vector of landmark observations
   */
  void DataAssociation(const std::vector<LandmarkObs> &predicted,
                       std::vector<LandmarkObs> &observations);

  /**
   * TransformToMapCoordinate to transform the sensor's measurements from its local sensor
   * coordinate system to the map's coordinate system
   * @param particle The given particle
   * @param observations Vector of landmark observations
   * @param transformed_obs Transformed vector of landmark observations
   */
  void TransformToMapCoordinate(const Particle &particle,
                                const std::vector<LandmarkObs> &observations,
                                std::vector<LandmarkObs> &transformed_obs);


  /**
   * FindInRangeLandmarks To find landmarks in sensor range for a given particle
   * @param sensor_range Sensor range
   * @param particle The given particle
   * @param map_landmarks Map class containing map landmarks
   * @param in_range_particles Vector of landmarks in sensor range for a given particle
   */
  void FindInRangeLandmarks(const double &sensor_range, const Particle &particle, const Map &map_landmarks,
                            std::vector<LandmarkObs> &in_range_particles);

  /**
   * CalculateParticleWeight To calculate the particle's weight
   * @param std_landmark[] Array of dimension 2 [standard deviation of x [m], standard deviation of y [m]]
   * @param predicted
   * @param observations Vector of landmark observations in MAP coordinate
   * @return A new value for particle's weight
   */
  double CalculateParticleWeight(const double *std_landmark,
                                 const std::vector<LandmarkObs> &predicted,
                                 const std::vector<LandmarkObs> &observations);

  /**
   * MultivariateGaussianProbability To calculate the Multivariate-Gaussian probability
   * @param x
   * @param y
   * @param mu_x
   * @param mu_y
   * @param std_x
   * @param std_y
   * @return
   */
  double MultivariateGaussianProbability(const double &x, const double &y, const double &mu_x,
                                         const double &mu_y, const double &std_x, const double &std_y);
  /**
   * UpdateWeights Updates the weights for each particle based on the likelihood of the
   * observed measurements.
   * @param sensor_range Range [m] of sensor
   * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
   * @param observations Vector of landmark observations
   * @param map Map class containing map landmarks
   */
  void UpdateWeights(const double &sensor_range, const double *std_landmark, const std::vector<LandmarkObs> &observations,
                     const Map &map_landmarks);

  /**
   * Resample Re-samples from the updated set of particles to form
   * the new set of particles.
   */
  void Resample();

  /**
   * Set a particles list of associations, along with the associations calculated world x,y coordinates
   * This can be a very useful debugging tool to make sure transformations are correct and associations
   * correctly connected
   */
  Particle SetAssociations(Particle &particle, const std::vector<int> &associations,
                           const std::vector<double> &sense_x, const std::vector<double> &sense_y);


  std::string GetAssociations(Particle best);

  std::string GetSenseX(Particle best);

  std::string GetSenseY(Particle best);

  /**
  * Initialized Returns whether particle filter is initialized yet or not.
  */
  const bool Initialized() const {
    return is_initialized;
  }
};


#endif /* PARTICLE_FILTER_H_ */
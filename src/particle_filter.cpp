/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <iostream> //rmiucic remove later
#include <string>//rmiucic remove later

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles

  std::default_random_engine gen;
  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std[0]);
  
  // This line Create normal distributions for y and theta
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i<num_particles; i++)
  {
    Particle tmp_particle;
    tmp_particle.x=dist_x(gen);
    tmp_particle.y=dist_y(gen);
    tmp_particle.theta=dist_theta(gen);
    tmp_particle.weight=1.0;
    particles.push_back(tmp_particle);
  }
  is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  std::default_random_engine gen;

  for (int i = 0; i<num_particles; i++)
  {
    if(yaw_rate==0.0)
    {
      particles[i].x=particles[i].x+velocity*delta_t*cos(particles[i].theta);
      particles[i].y=particles[i].y+velocity*delta_t*sin(particles[i].theta);
      particles[i].theta=particles[i].theta;
    }
    else
    {//move each particle filter according to bicycle model
      particles[i].x=particles[i].x+velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
      particles[i].y=particles[i].y+velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
      particles[i].theta=particles[i].theta+yaw_rate*delta_t;
    }
    /*Add Gaussian Noise around mean of the calculated particle position */
    // This line creates a normal (Gaussian) distribution for x
    normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    // This line Create normal distributions for y and theta
    normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
    //particle filter x,y, and theta is draw from gausian distribution
    particles[i].x=dist_x(gen);
    particles[i].y=dist_y(gen);
    particles[i].theta=dist_theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  /* predicted - prediction measurements between one particular particle and all 
   *             of the map landmarks within sensor range 
   * observations - actual landmark measurements gathered from the LIDAR
   * 
   * This function performs nearest neighbor data associations and assign each 
   * sensor observation the map landmark ID associated with it
   */
    for(unsigned int l=0;l<observations.size();l++)
    {
      double min_dist=100.0; //some large number
      for(unsigned int m=0;m<predicted.size();m++)
      {
        double dist_obs_prd=dist(observations[l].x,observations[l].y,predicted[m].x,predicted[m].y);
        if (dist_obs_prd<min_dist)
        {
          min_dist=dist_obs_prd;
          observations[l].id=predicted[m].id; //update observation
        }
      }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  for ( int i = 0; i<num_particles; i++)
  {
    //Convert observations to world coordinate system
    LandmarkObs obs_map; //temporary variable
    vector<LandmarkObs> observations_map;
    for(unsigned int j=0;j<observations.size();j++)
    {
      obs_map.x=particles[i].x+(cos(particles[i].theta)*observations[j].x)
                 -(sin(particles[i].theta) * observations[j].y);
      obs_map.y=particles[i].y+(sin(particles[i].theta)*observations[j].x)
                 +(cos(particles[i].theta) * observations[j].y);

      observations_map.push_back(obs_map);
    }

    //select predicted landmarks that are in range of the sensor 
    vector<LandmarkObs> predicted_map;
    LandmarkObs prd_map; //temporary variable
    for(unsigned int k=0;k<map_landmarks.landmark_list.size();k++)
    {
      prd_map.id=map_landmarks.landmark_list[k].id_i;
      prd_map.x=map_landmarks.landmark_list[k].x_f;
      prd_map.y=map_landmarks.landmark_list[k].y_f;
      if( dist(prd_map.x,prd_map.y, particles[i].x, particles[i].y) < sensor_range)
      {
        predicted_map.push_back(prd_map);
      }
    }
    //Associate each observation to coresponding landmark
    dataAssociation(predicted_map,observations_map);
    vector<int> associations_tmp;
    vector<double> sense_x_tmp;
    vector<double> sense_y_tmp;
    vector<double> obs_prob;
    for(unsigned int l=0;l<observations_map.size();l++)
    {
      associations_tmp.push_back(observations_map[l].id);
      sense_x_tmp.push_back(observations_map[l].x);
      sense_y_tmp.push_back(observations_map[l].y);
      
      for (unsigned int m=0; m<predicted_map.size(); m++)
      {
        if(predicted_map[m].id==observations_map[l].id)
        {        
          obs_prob.push_back(multiv_prob(std_landmark[0],
                                  std_landmark[1], 
                                  observations_map[l].x, 
                                  observations_map[l].y,
                                  predicted_map[m].x,
                                  predicted_map[m].y));
          break;
        }

      }
    }
    double particle_prob=1;
    for(unsigned int n=0;n<obs_prob.size();n++)
    {
      particle_prob=obs_prob[n];
    }
    particles[i].weight=particle_prob;
    SetAssociations(particles[i],associations_tmp, sense_x_tmp,sense_y_tmp);
  }//end(int i = 0; i<num_particles; i++)
  //normalize particle weights 
  double weights_sum=0;
  for(unsigned int p=0;p<particles.size();p++)
  {
    weights_sum=weights_sum+particles[p].weight;
  }
  for(unsigned int p=0;p<particles.size();p++)
  {
    particles[p].weight=particles[p].weight/weights_sum;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<double> weights_all;
  for(unsigned int p=0;p<particles.size();p++)
  {
    weights_all.push_back(particles[p].weight);
  }

  std::default_random_engine generator;
  std::discrete_distribution<int> distribution (weights_all.begin(),weights_all.end());
  
  vector<Particle> particles_tmp;
  for(unsigned int p=0;p<particles.size();p++)
  {
    particles_tmp.push_back(particles[distribution(generator)]);
  }
  particles=particles_tmp;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

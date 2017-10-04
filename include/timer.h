// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_TIMER_H_INCLUDED
#define MPS_TIMER_H_INCLUDED

#include <chrono>
#include <ctime>
#include <iostream>
#include <boost/format.hpp>

namespace tiny_mps {

// Holds data on time
class Timer {
 public:
  explicit Timer(const Condition & condition) {
    initialize(condition);
  }
  // Timer is neither copyable nor movable.
  Timer(const Timer&) = delete;
  Timer& operator=(const Timer&) = delete;
  virtual ~Timer(){}

  // Initializes data members using Condition reference.
  // Stores the current time point.
  inline void initialize(const Condition& condition) {
    this->current_time = condition.initial_time;
    this->initial_time = condition.initial_time;
    this->finish_time = condition.finish_time;
    this->current_delta_time = condition.delta_time;
    this->initial_delta_time = condition.delta_time;
    this->min_delta_time = condition.min_delta_time * condition.delta_time;
    this->output_interval = condition.output_interval;
    this->next_output_time = condition.initial_time;
    this->loop_count = 0;
    this->output_count = 0;
    start_chrono = std::chrono::system_clock::now();
    std::time_t start = std::chrono::system_clock::to_time_t(start_chrono);
    std::cout << "Started timer at " << std::ctime(&start) << std::endl;
  }

  inline bool hasNextLoop() const {
    if (current_time > finish_time) return false;
    else return true;
  }

  inline void update() {
    if (isOutputTime()) {
        next_output_time += output_interval;
        ++output_count;
    }
    current_time += current_delta_time;
    ++loop_count;
  }

  inline void printCompuationTime() {
    using std::chrono::duration_cast;
    using std::chrono::system_clock;
    using std::chrono::hours;
    using std::chrono::minutes;
    using std::chrono::seconds;
    using std::chrono::milliseconds;
    auto end = system_clock::now();
    auto dur = end - start_chrono;
    std::cout << boost::format("Computation Time: %02dh %02dmin %02d.%03ds.")
        % duration_cast<hours>(dur).count()
        % (duration_cast<minutes>(dur).count() % 60)
        % (duration_cast<seconds>(dur).count() % 60)
        % (duration_cast<milliseconds>(dur).count() % 1000) << std::endl;
  }

  inline void printTimeInfo() {
    std::cout << boost::format("Time step: %06d, Current time: %e, Delta time: %e")
        % getLoopCount()
        % getCurrentTime()
        % getCurrentDeltaTime() << std::endl;
  }

  inline void limitCurrentDeltaTime(double max_speed, const Condition& condition) {
    if (max_speed <= 0) return;
    current_delta_time = initial_delta_time;
    double dt = condition.average_distance * condition.courant_number / max_speed;
    current_delta_time = std::min(dt, current_delta_time);
    if (condition.viscosity_calculation == false) return;
    dt = condition.diffusion_number * condition.average_distance * condition.average_distance
        / condition.kinematic_viscosity;
    current_delta_time = std::min(dt, current_delta_time);
  }

  inline bool isUnderMinDeltaTime() {
    return getCurrentDeltaTime() < min_delta_time;
  }
  inline bool isOutputTime() const {
    return current_time >= next_output_time;
  }

  inline double getCurrentTime() const { return current_time; }
  inline double getInitialTime() const { return initial_time; }
  inline double getFinishTime() const { return finish_time; }
  inline double getCurrentDeltaTime() const { return current_delta_time; }
  inline double getInitialDeltaTime() const { return initial_delta_time; }
  inline void setInitialDeltaTime(double delta_time) { initial_delta_time = delta_time; }
  inline int getLoopCount() const { return loop_count; }
  inline int getOutputCount() const { return output_count; }

private:
  // Stores the time point when called initialize.
  std::chrono::system_clock::time_point start_chrono;

  double current_time;
  double initial_time;
  double finish_time;
  double current_delta_time;
  double initial_delta_time;
  double min_delta_time;
  double output_interval;
  // The next time point when to write output files.
  double next_output_time;
  int loop_count;
  int output_count;
};

} // namespace tiny_mps
#endif //MPS_TIMER_H_INCLUDED
// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_TIMER_H_INCLUDED
#define MPS_TIMER_H_INCLUDED

#include <boost/format.hpp>

namespace tiny_mps {

// Holds data on time
class Timer {
public:
    Timer(const Condition & condition) {
        initialize(condition);
    }

    virtual ~Timer(){}

    inline void initialize(const Condition& condition) {
        this->current_time = condition.initial_time;
        this->initial_time = condition.initial_time;
        this->finish_time = condition.finish_time;
        this->current_delta_time = condition.delta_time;
        this->initial_delta_time = condition.delta_time;
        this->min_delta_time = condition.min_delta_time * condition.delta_time;
        this->output_interval = condition.output_interval;
        this->loop_count = 0;
        this->output_count = 0;
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
    
    inline void printTimeInfo() {
        std::cout << boost::format("Time step: %08d, Current time: %f, Delta time: %f")
        % getLoopCount() % getCurrentTime() % getCurrentDeltaTime() << std::endl;
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
    double current_time;
    double initial_time;
    double finish_time;
    double current_delta_time;
    double initial_delta_time;
    double min_delta_time;
    double output_interval;
    double next_output_time;
    int loop_count;
    int output_count;
};

} // namespace tiny_mps
#endif //MPS_TIMER_H_INCLUDED
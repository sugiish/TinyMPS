#ifndef MPS_TIMER_H_INCLUDED
#define MPS_TIMER_H_INCLUDED

class Timer
{
public:
    Timer();
    virtual ~Timer();

    inline void initialize()
    {
        current_time = 0.0;
        initial_time = 0.0;
        finish_time = 0.0;

        current_delta_time = 0.0;
        initial_delta_time = 0.0;

        cfl_condition = 0.2;
    }

    inline void update()
    {
        current_time += current_delta_time;
    }

    inline void limitCurrentDeltaTime(double max_speed, double particle_distance, double diffusion_number, double kinematic_viscosity)
    {
        current_delta_time = initial_delta_time;
        
        if (max_speed <= 0) return;

        double dt = particle_distance * cfl_condition / max_speed;
        if (dt < current_delta_time)
        {
            current_delta_time = dt;
        }

        dt = diffusion_number * particle_distance * particle_distance / kinematic_viscosity;
        if (dt < current_delta_time)
        {
            current_delta_time = dt;
        }
    }

    inline double getCurrentTime()
    {
        return current_time;
    }

    inline double getInitialTime()
    {
        return initial_time;
    }

    inline double getFinishTime()
    {
        return finish_time;
    }

    inline double getCurrentDeltaTime()
    {
        return current_delta_time;
    }

    inline double getInitialDeltaTime()
    {
        return initial_delta_time;
    }

    inline void setInitialDeltaTime(double delta_time)
    {
        initial_delta_time = delta_time;
    }

    inline double getCFLCondition()
    {
        return cfl_condition;
    }

    inline void setCFLCondition(double value)
    {
        cfl_condition = value;
    }

private:
    double current_time;
    double initial_time;
    double finish_time;

    double current_delta_time;
    double initial_delta_time;

    double cfl_condition;

};

#endif //MPS_TIMER_H_INCLUDED
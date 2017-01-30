#include <functional>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <cmath>
#include <cassert>


// Simulation constants
static const float              MAX_SIMULATION_LENGTH = 7200.0f;
static float                    ARRIVIAL_RATE = 1.0f/10.0f;
static float                    UNION_REGULATION = 12.0f;
static std::pair<float, float>  UNLOAD_WORK_RANGE(3.5f, 4.5f);
static std::pair<float, float>  CREW_REMAINING_TIME_RANGE(6.0f, 11.0f);
static std::pair<float, float>  CREW_REPLACEMENT_TIME_RANGE(2.5f, 3.5f);


// Utils.
class UniformDist
{
public:
        UniformDist();
        UniformDist(uint32_t seed);
        float   next_dice(float a, float b) const;
};

UniformDist::UniformDist()
{
}

UniformDist::UniformDist(uint32_t seed)
{
        std::srand(seed);
}

float UniformDist::next_dice(float a, float b) const
{
        return a + ((float) std::rand()/RAND_MAX)*(b - a);
}

class PoissonProcess
{
public:
        PoissonProcess();
        PoissonProcess(uint32_t seed);
        float   next_arrival(float rate) const;
private:
        UniformDist     m_u;
};

PoissonProcess::PoissonProcess()
{
}

PoissonProcess::PoissonProcess(uint32_t seed):
        m_u(seed)
{
}

float PoissonProcess::next_arrival(float r) const
{
        return -1.0f/r*std::log(m_u.next_dice(0, 1));
}

// Entity storage.
class Entity {};

class UnloadingSystem: public Entity
{
public:
        uint32_t        spawn_new_train(float t, float r, float& next_spawn, float a, float b, float& work_amount);
        float&          work_amount(uint32_t tid);
        float           kill_train(uint32_t tid, float t);

        uint32_t        train_spawned() const;
        float           avg_tis() const;
        float           max_tis() const;
private:
        struct TrainStats
        {
                float   ts;
                float   work_amount;

                TrainStats() {};
                TrainStats(float ts, float work_amount):
                        ts(ts), work_amount(work_amount) {}
        };

        std::map<uint32_t, TrainStats>  m_train_in;
        uint32_t                        m_num_trains = 0;
        float                           m_max_tis = 0.0f;
        float                           m_sum_tis = 0.0f;
        PoissonProcess                  m_p;
        UniformDist                     m_u;
};

uint32_t UnloadingSystem::UnloadingSystem::spawn_new_train(float t, float r, float& next_spawn,
                                                           float a, float b, float &work_amount)
{
        next_spawn = m_p.next_arrival(r);
        work_amount = m_u.next_dice(a, b);

        m_num_trains ++;
        m_train_in.insert(std::pair<uint32_t, TrainStats>(m_num_trains, TrainStats(t, work_amount)));
        return m_num_trains;
}

float& UnloadingSystem::work_amount(uint32_t tid)
{
        return m_train_in[tid].work_amount;
}

float UnloadingSystem::kill_train(uint32_t tid, float t)
{
        float tis = t - m_train_in[tid].ts;
        m_max_tis = std::max(m_max_tis, tis);
        m_sum_tis += tis;
        m_train_in.erase(tid);
        return tis;
}

uint32_t UnloadingSystem::train_spawned() const
{
        return m_num_trains;
}

float UnloadingSystem::avg_tis() const
{
        return m_sum_tis/m_num_trains;
}

float UnloadingSystem::max_tis() const
{
        return m_max_tis;
}

class TrainQueue: public Entity
{
public:
        void                            train_arrive(uint32_t tid, float t);
        bool                            unload_next_train(float t);
        bool                            peek_next_train(uint32_t& tid);
        const std::list<uint32_t>&      all_trains() const;

        uint32_t        max_length() const;
        float           time_average() const;
private:
        void            update(float t);
        std::list<uint32_t>     m_trains;
        uint32_t                m_max_length = 0;
        float                   m_moment = 0.0f;
        float                   m_t = 0.0f;
};

void TrainQueue::update(float t)
{
        assert(t >= m_t);
        m_moment += (t - m_t)*m_trains.size();
        m_t = t;
}

void TrainQueue::train_arrive(uint32_t tid, float t)
{
        update(t);

        m_trains.push_back(tid);
        m_max_length = std::max(m_max_length, (uint32_t) m_trains.size());
}

bool TrainQueue::unload_next_train(float t)
{
        if (m_trains.empty())
                return false;
        update(t);

        m_trains.pop_front();
        return true;
}

bool TrainQueue::peek_next_train(uint32_t& tid)
{
        if (m_trains.empty())
                return false;
        tid = m_trains.front();
        return true;
}

const std::list<uint32_t>& TrainQueue::all_trains() const
{
        return m_trains;
}

uint32_t TrainQueue::max_length() const
{
        return m_max_length;
}

float TrainQueue::time_average() const
{
        if (m_t == 0)
                return 0;
        else
                return m_moment/m_t;
}

class Dock: public Entity
{
public:
        void            unload(uint32_t tid, float t);
        void            hog_out(float t);
        void            wait(float t);

        uint32_t        current_train() const;
        bool            is_busy() const;
        bool            is_hogged_out() const;
        bool            is_waiting() const;

        float           p_busy() const;
        float           p_idle() const;
        float           p_hogged_out() const;
private:
        uint32_t        m_status = 0;
        uint32_t        m_curr_train = 0;
        float           m_dock_time[3] = {0.0f};
        float           m_t = 0.0f;
};

void Dock::unload(uint32_t tid, float t)
{
        m_dock_time[m_status] += (t - m_t);
        m_t = t;
        m_curr_train = tid;
        m_status = 1;
}

void Dock::hog_out(float t)
{
        m_dock_time[m_status] += (t - m_t);
        m_t = t;
        m_status = 2;
}

void Dock::wait(float t)
{
        m_dock_time[m_status] += (t - m_t);
        m_t = t;
        m_status = 0;
}

uint32_t Dock::current_train() const
{
        return m_curr_train;
}

bool Dock::is_busy() const
{
        return m_status == 1;
}

bool Dock::is_hogged_out() const
{
        return m_status == 2;
}

bool Dock::is_waiting() const
{
        return m_status == 0;
}

float Dock::p_busy() const
{
        return m_dock_time[1]/(m_dock_time[0] + m_dock_time[1] + m_dock_time[2]);
}

float Dock::p_idle() const
{
        return m_dock_time[0]/(m_dock_time[0] + m_dock_time[1] + m_dock_time[2]);
}

float Dock::p_hogged_out() const
{
        return m_dock_time[2]/(m_dock_time[0] + m_dock_time[1] + m_dock_time[2]);
}

class Crews: public Entity
{
public:
        float                   spawn_crew(uint32_t tid, float a, float b);
        bool                    work(uint32_t tid, float amount, float& t_worked);
        float                   get_arrival_time(uint32_t tid, float t) const;
        float                   replace(uint32_t tid, float t, float a, float b);
        void                    kill_crew(uint32_t tid);

        void                    hist(std::map<uint32_t, uint32_t>& hist) const;
private:
        struct CrewStats
        {
                float   t_remain;
                float   t_arrive;

                CrewStats() {}
                CrewStats(float t_remain, float t_arrive):
                        t_remain(t_remain), t_arrive(t_arrive) {}
        };

        std::map<uint32_t, CrewStats>   m_crews;
        std::map<uint32_t, uint32_t>    m_num_hogged_out;
        UniformDist                     m_u;
};

float Crews::spawn_crew(uint32_t tid, float a, float b)
{
        float remain = m_u.next_dice(a, b);
        m_crews.insert(std::pair<uint32_t, CrewStats>(tid, CrewStats(remain, -1.0f)));
        m_num_hogged_out.insert(std::pair<uint32_t, uint32_t>(tid, 0));
        return remain;
}

bool Crews::work(uint32_t tid, float amount, float& t_worked)
{
        if (m_crews[tid].t_remain == 0) {
                // Unable to work when hogged-out.
                t_worked = 0.0f;
                return false;
        }
        float remain = m_crews[tid].t_remain - amount;
        if (remain < 0) {
                // Hogged out.
                t_worked = amount + remain;
                m_crews[tid].t_remain = 0;
                m_num_hogged_out[tid] ++;
                return false;
        } else {
                t_worked = amount;
                m_crews[tid].t_remain = remain;
                return true;
        }
}

float Crews::get_arrival_time(uint32_t tid, float t) const
{
        if (t > m_crews.at(tid).t_arrive)
                return -1;
        else
                return m_crews.at(tid).t_arrive;
}

float Crews::replace(uint32_t tid, float t, float a, float b)
{
        assert(m_crews[tid].t_remain == 0);
        float travel_time = m_u.next_dice(a, b);
        m_crews[tid].t_arrive = t + travel_time;
        m_crews[tid].t_remain = ::UNION_REGULATION - travel_time;
        return travel_time;
}

void Crews::kill_crew(uint32_t tid)
{
        m_crews.erase(tid);
}

void Crews::hist(std::map<uint32_t, uint32_t>& hist) const
{
        for (std::map<uint32_t, uint32_t>::const_iterator it = m_num_hogged_out.begin();
             it != m_num_hogged_out.end(); ++ it) {
                uint32_t n_hogged = it->second;
                std::map<uint32_t, uint32_t>::iterator h_it = hist.find(n_hogged);
                if (h_it == hist.end()) {
                        hist.insert(std::pair<uint32_t, uint32_t>(n_hogged, 1));
                } else {
                        h_it->second ++;
                }
        }
}


// Entities.
class Entities
{
public:
        UnloadingSystem         unloading_system;
        TrainQueue              train_queue;
        Dock                    dock;
        Crews                   crews;
public:
};

// Events.
class EventQueue;
class IEvent
{
public:
        IEvent(float when);
        virtual void hit(EventQueue& e, Entities& entities, float t, float tl) = 0;
        bool operator<(const IEvent& e) const;
        bool operator>(const IEvent& e) const;
public:
        float   when;
};

IEvent::IEvent(float when)
{
        this->when = when;
}

bool IEvent::operator<(const IEvent& e) const
{
        return when < e.when;
}

bool IEvent::operator>(const IEvent& e) const
{
        return when > e.when;
}

class Spawn: public IEvent
{
public:
        Spawn(float when): IEvent(when) {}
        void hit(EventQueue& e, Entities& entities, float t, float tl) override;
};

class Kill: public IEvent
{
public:
        Kill(float when, uint32_t which): IEvent(when), m_which(which) {}
        void hit(EventQueue& e, Entities& entities, float t, float tl) override;
private:
        uint32_t        m_which;
};

class Unload: public IEvent
{
public:
        Unload(float when, uint32_t which):
                IEvent(when), m_which(which) {}
        void hit(EventQueue& e, Entities& entities, float t, float tl) override;
private:
        uint32_t        m_which;
};

class QueueWait: public IEvent
{
public:
        QueueWait(float when):
                IEvent(when) {}
        void hit(EventQueue& e, Entities& entities, float t, float tl) override;
};

class Swap: public IEvent
{
public:
        Swap(float when, float leftover, uint32_t which):
                IEvent(when), m_leftover(leftover), m_which(which) {}
        void hit(EventQueue& e, Entities& entities, float t, float tl) override;
private:
        float           m_leftover;
        uint32_t        m_which;
};

class Overtime: public IEvent
{
public:
        Overtime(float when, uint32_t which):
                IEvent(when), m_which(which) {}
        void hit(EventQueue& e, Entities& entities, float t, float tl) override;
private:
        uint32_t        m_which;
};

class Replacement: public IEvent
{
public:
        Replacement(float when, uint32_t which):
                IEvent(when), m_which(which) {}
        void hit(EventQueue& e, Entities& entities, float t, float tl) override;
private:
        uint32_t        m_which;
};

// Event queue.
class EventQueue
{
public:
        void    schedule(IEvent *e);
        IEvent* pop();
        bool    has_done() const;

private:
        class EventComparator
        {
        public:
                bool operator()(const IEvent* e0, const IEvent* e1) const
                {
                        return *e0 > *e1;
                }
        };
        std::priority_queue<IEvent*,
                            std::vector<IEvent*>,
                            EventComparator>    m_events;
};

void EventQueue::schedule(IEvent* e)
{
        m_events.push(e);
}

IEvent* EventQueue::pop()
{
        IEvent* e = m_events.top();
        m_events.pop();
        return e;
}

bool EventQueue::has_done() const
{
        return m_events.empty();
}

// Event defns.
void Spawn::hit(EventQueue& e, Entities& entities, float t, float tl)
{
        float next_spawn, work_amount;
        uint32_t tid = entities.unloading_system.spawn_new_train(t, ::ARRIVIAL_RATE, next_spawn,
                                ::UNLOAD_WORK_RANGE.first, ::UNLOAD_WORK_RANGE.second, work_amount);
        if (tl > 0.0f) {
                // Next arrvial events should occur.
                e.schedule(new Spawn(t + next_spawn));
        } else {
                // Stop source and finish the remaining events.
        }
        float remaining = entities.crews.spawn_crew(tid,
                                                    ::CREW_REMAINING_TIME_RANGE.first,
                                                    ::CREW_REMAINING_TIME_RANGE.second);
        entities.train_queue.train_arrive(tid, t);
        e.schedule(new QueueWait(t + 0));
        std::cout << "At t=" << t << " hrs -> ";
        std::cout << "Train " << tid
                  << " arrives with anticipated work amount of " << work_amount << " hours and"
                  << " crews' remaing hours of " << remaining
                  << std::endl;
}


void Kill::hit(EventQueue& e, Entities& entities, float t, float tl)
{
        float tis = entities.unloading_system.kill_train(m_which, t);
        entities.crews.kill_crew(m_which);
        entities.dock.wait(t);
        e.schedule(new QueueWait(t + 0));

        std::cout << "At t=" << t << " hrs -> ";
        std::cout << "Train " << m_which
                  << " departs with spending " << tis << " hours in system."
                  << std::endl;
}

void QueueWait::hit(EventQueue& e, Entities& entities, float t, float tl)
{
        if (entities.dock.is_waiting()) {
                // Unload next train.
                uint32_t next_train;
                if (!entities.train_queue.peek_next_train(next_train))
                        // No trains in queue.
                        return;
                else {
                        float a = entities.crews.get_arrival_time(next_train, t);
                        if (a > t) {
                                std::cout << "At t=" << t << " hrs -> ";
                                std::cout << "Train " << next_train
                                          << " is blocking in queue while there are "
                                          << entities.train_queue.all_trains().size() << " trains in queue."
                                          << std::endl;
                        } else {
                                entities.train_queue.unload_next_train(t);
                                e.schedule(new Unload(t + 0, next_train));
                                std::cout << "At t=" << t << " hrs -> ";
                                std::cout << "Train " << next_train
                                          << " will leave the queue while there are "
                                          << entities.train_queue.all_trains().size() << " trains in queue."
                                          << std::endl;
                        }
                }
        } else {
                std::cout << "At t=" << t << " hrs -> ";
                std::cout << "There are "
                          << entities.train_queue.all_trains().size() << " trains in queue."
                          << std::endl;
        }
}

void Unload::hit(EventQueue& e, Entities& entities, float t, float tl)
{
        // Starting unloading.
        float& total = entities.unloading_system.work_amount(m_which);
        float worked;
        if (!entities.crews.work(m_which, total, worked)) {
                // Hogged out.
                e.schedule(new Overtime(t + worked, m_which));
        } else {
                // Complete.
                e.schedule(new Kill(t + worked, m_which));
        }
        entities.dock.unload(m_which, t);
        total -= worked;

        std::cout << "At t=" << t << " hrs -> ";
        std::cout << "Unloading train " << m_which
                  << " with the work amount left " << total << " hours."
                  << std::endl;

        // Trains in queue are also working.
        for (uint32_t tid: entities.train_queue.all_trains()) {
                float w = 0.0f, leftover;
                if (!entities.crews.work(tid, worked, w)) {
                        e.schedule(new Swap(t + w, leftover, tid));
                }
        }
}

void Swap::hit(EventQueue& e, Entities& entities, float t, float tl)
{
        std::cout << "At t=" << t << " hrs -> ";
        std::cout << "Train " << m_which
                  << " waiting in queue needs to swap crew members."
                  << std::endl;
        float w;
        if (!entities.crews.work(m_which, m_leftover, w)) {
                float a = entities.crews.get_arrival_time(m_which, t);
                if (a < 0) {
                        a = t + entities.crews.replace(m_which, t,
                                                       ::CREW_REPLACEMENT_TIME_RANGE.first,
                                                       ::CREW_REPLACEMENT_TIME_RANGE.second);
                }
                float leftover = m_leftover - (a - t);
                if (leftover > 0)
                        e.schedule(new Swap(a, leftover, m_which));
                else
                        e.schedule(new QueueWait(a));
        } else
                e.schedule(new QueueWait(t));
}

void Overtime::hit(EventQueue& e, Entities& entities, float t, float tl)
{
        entities.dock.hog_out(t);
        float travel_time = entities.crews.replace(m_which, t,
                                                   ::CREW_REPLACEMENT_TIME_RANGE.first,
                                                   ::CREW_REPLACEMENT_TIME_RANGE.second);
        e.schedule(new Replacement(t + travel_time, m_which));
        std::cout << "At t=" << t << " hrs -> ";
        std::cout << "Crews of train " << m_which
                  << " needs to be replaced."
                  << std::endl;
}


void Replacement::hit(EventQueue& e, Entities& entities, float t, float tl)
{
        e.schedule(new Unload(t + 0, m_which));
        std::cout << "At t=" << t << " hrs -> ";
        std::cout << "Crews of train " << m_which
                  << " has been replaced."
                  << std::endl;
}


// Simulation.
class Simulation
{
public:
        Simulation(const std::vector<IEvent*>& initial_events,
                   Entities& entities);
        void    simulate(float length);
private:

        Entities&       m_entities;
        EventQueue      m_events;
};

Simulation::Simulation(const std::vector<IEvent*>& initial_events,
                       Entities& entities):
        m_entities(entities)
{
        for (IEvent* event: initial_events)
                m_events.schedule(event);
}

void Simulation::simulate(float length)
{
        while (!m_events.has_done()) {
                IEvent* e = m_events.pop();
                e->hit(m_events, m_entities, e->when, length - e->when);
                delete e;
        }
}

int main(int argc, char *argv[])
{
        std::srand(time(0));

        std::vector<IEvent*> initial_events;
        initial_events.push_back(new Spawn(0.0f));

        Entities entities;
        Simulation s(initial_events, entities);

        std::cout << "=================================================" << std::endl
                  << "Simulation begins " << std::endl
                  << "=================================================" << std::endl;
        s.simulate(::MAX_SIMULATION_LENGTH);

        std::cout << std::endl;
        std::cout << "=================================================" << std::endl
                  << "Results: " << std::endl
                  << "=================================================" << std::endl;

        std::cout << "Trains" << std::endl;
        std::cout << entities.unloading_system.train_spawned() << " trains served." << std::endl;
        std::cout << "Average time-in-system: " << entities.unloading_system.avg_tis() << " hours." << std::endl;
        std::cout << "Peak time-in-system: " << entities.unloading_system.max_tis() << " hours." << std::endl;

        std::cout << std::endl;
        std::cout << "Dock" << std::endl;
        std::cout << "Busy: " << 100.0f*entities.dock.p_busy() << "%" << std::endl;
        std::cout << "Hogged-out: " << 100.0f*entities.dock.p_hogged_out() << "%" << std::endl;
        std::cout << "Idle: " << 100.0f*entities.dock.p_idle() << "%" << std::endl;

        std::cout << std::endl;
        std::cout << "Queueing" << std::endl;
        std::cout << "Time average: " << entities.train_queue.time_average() << " trains in queue." << std::endl;
        std::cout << "Peak: " << entities.train_queue.max_length() << " trains in queue." << std::endl;

        std::cout << std::endl;
        std::cout << "Hogged-out" << std::endl;
        std::map<uint32_t, uint32_t> hog_hist;
        entities.crews.hist(hog_hist);
        uint32_t max_hog = entities.unloading_system.train_spawned();

        for (uint32_t i = 0; i < hog_hist.size(); i ++) {
                std::cout << i << " times: ";

                std::map<uint32_t, uint32_t>::iterator h_it = hog_hist.find(i);
                if (h_it != hog_hist.end()) {
                        uint32_t h = h_it->second;
                        float r = (float) h/max_hog;
                        uint32_t n = (uint32_t) (r*100.0f);

                        for (uint32_t j = 0; j < n; j ++)
                                std::cout << '+';

                        std::cout << "  " << h << std::endl;
                } else {
                        std::cout << "  " << 0 << std::endl;
                }
        }
        return 0;
}

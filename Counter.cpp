#include "Counter.h"

int Counter::count = 0;

Counter::Counter() 
    : state(CounterState::Inactive)
{
    ++count;
}

int Counter::get_count() {
    return count;
}

void Counter::start() {
    state = CounterState::Active;
}

void Counter::finish() {
    state = CounterState::Finished;
}

CounterState Counter::get_state() const {
    return state;
}
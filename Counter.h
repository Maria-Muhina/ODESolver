#pragma once

enum class CounterState {
    Inactive,
    Active,
    Finished
};

class Counter {
    static int count;
    CounterState state;

    public:
        Counter();
        static int get_count();
        void start();
        void finish();
        CounterState get_state() const;
};

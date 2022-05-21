#include "nbody.h"

#include <iostream>

int main(int argc, char ** argv)
{
    if (argc == 4) {
        try {
            std::size_t time = std::stoull(std::string(argv[2]));
            std::size_t time_step = std::stoull(std::string(argv[3]));

            FastPositionTracker tracker = FastPositionTracker(argv[1]);
            tracker.track(tracker.bodies().front().get_m_name(), time, time_step);
            for (const Body & i : tracker.bodies()) {
                std::cout << i.get_m_name() << '\t' << i.get_m_pos().x << '\t' << i.get_m_pos().y << std::endl;
            }
            return 0;
        }
        catch (std::logic_error const & ex) {
            std::cerr << "Argument parse error: " << ex.what() << std::endl;
        }
    }
    std::cerr << "Syntax: \nnbody BODY_FILE END_TIME TIME_STEP" << std::endl;
    return 1;
}

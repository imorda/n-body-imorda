#include <iostream>
#include <string_view>

struct Cartesian
{
    double x;
    double y;
};

// Single body representation, required for Problem 1 and Problem 2
class Body
{
public:
    Body(...);

    double distance(const Body & b);

    // calculate the force-on current body by the 'b' and add the value to accumulated force value
    void add_force(const Body & b);
    // reset accumulated force value
    void reset_force();

    // update body's velocity and position
    void update(double delta_t);

    friend std::ostream & operator<<(std::ostream &, const Body &);

    // The methods below to be done for Burnes-Hut algorithm only
    // Test if body is in quadrant
    bool in(const Quadrant q);
    // Create new body representing center-of-mass of the invoking body and 'b'
    Body plus(const Body & b);
};

// Quadrant representation, required for Problem 2
class Quadrant
{
public:
    // Create quadrant with center (x, y) and size 'lenth'
    Quadrant(Cartesian center, double length);

    // Test if point (x, y) is in the quadrant
    bool contains(Cartesian p);
    double length();

    // The four methods below construct new Quadrant representing sub-quadrant of the invoking quadrant
    Quadrant nw();
    Quadrant ne();
    Quadrant sw();
    Quadrant se();

    friend std::ostream & operator<<(std::ostream &, const Quadrant &);
};

// Burnes-Hut tree representation, required for Problem 2
class BHTreeNode
{
public:
    BHTreeNode(...);

    void insert(Body b);
    // Update net acting force-on 'b'
    void update_force(Body & b);
};

class PositionTracker
{
protected:
    PositionTracker(const std::string & filename);

public:
    virtual Track track(const std::string & body_name, size_t end_time, size_t time_step) = 0;
};

class BasicPositionTracker : public PositionTracker
{
public:
    BasicPositionTracker(const std::string & filename);
    Track track(const std::string & body_name, size_t end_time, size_t time_step) override;
};

class FastPositionTracker : public PositionTracker
{
public:
    FastPositionTracker(const std::string & filename);
    Track track(const std::string & body_name, size_t end_time, size_t time_step) override;
};

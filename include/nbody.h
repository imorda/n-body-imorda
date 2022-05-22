#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <string_view>
#include <vector>

namespace constants {
const double G = (6.67e-11);
const double THETA = 0.5;
} // Namespace constants

struct Cartesian
{
    double x;
    double y;
};

// Quadrant representation, required for Problem 2
class Quadrant
{
private:
    Cartesian m_center;
    double m_length;

public:
    // Create quadrant with center (x, y) and size 'length'
    Quadrant(Cartesian center, double length);
    // Test if point (x, y) is in the quadrant
    bool contains(const Cartesian & p) const;

    double length() const;
    // The four methods below construct new Quadrant representing sub-quadrant of the invoking quadrant
    Quadrant nw() const;

    Quadrant ne() const;
    Quadrant sw() const;
    Quadrant se() const;

    friend std::ostream & operator<<(std::ostream &, const Quadrant &);
};

// Single body representation, required for Problem 1 and Problem 2
class Body
{
private:
    Cartesian m_pos;
    Cartesian m_speed;
    Cartesian m_force = {0, 0};
    double m_mass;
    std::string m_name;

public:
    Body(const Cartesian & m_pos, const Cartesian & m_speed, double m_mass, const std::string & mName);
    Body();
    double distance(const Body & b) const;
    const Cartesian & get_m_pos() const;
    const std::string & get_m_name() const;
    // calculate the force-on current body by the 'b' and add the value to accumulated force value
    void add_force(const Body & b);
    // reset accumulated force value
    void reset_force();

    // update body's velocity and position
    void update(double delta_t);

    friend std::ostream & operator<<(std::ostream &, const Body &);

    // The methods below to be done for Burnes-Hut algorithm only
    // Test if body is in quadrant
    bool in(const Quadrant & q) const;
    // Create new body representing center-of-m_mass of the invoking body and 'b'
    Body plus(const Body & b) const;
};

using Track = std::vector<Cartesian>;

class PositionTracker
{
protected:
    double m_size;
    std::vector<Body> m_bodies;
    PositionTracker(const std::string & filename);
    virtual void track_impl() = 0;

public:
    Track track(const std::string & body_name, size_t end_time, size_t time_step);
    const std::vector<Body> & bodies();
    virtual ~PositionTracker() = default;
};

class BasicPositionTracker : public PositionTracker
{
protected:
    void track_impl() override;

public:
    BasicPositionTracker(const std::string & filename);
};

// Burnes-Hut tree representation, required for Problem 2
class BHTreeNode
{
protected:
    Body m_representative;

public:
    BHTreeNode(const Body & b);
    virtual ~BHTreeNode() = default;
    const Body & get_m_representative() const;
    virtual void insert(const Body &);
    // Update net acting force-on 'b'
    virtual void update_force(Body & b);

protected:
    BHTreeNode() = default;
    virtual bool is_far_enough(const Body & b);
};
class BHTreeInternal : public BHTreeNode
{
private:
    std::array<std::array<std::unique_ptr<BHTreeNode>, 2>, 2> m_children;
    Quadrant m_area;

public:
    BHTreeInternal(const Quadrant & mArea);
    void insert(const Body & b) override;
    void update_force(Body & b) override;

protected:
    bool is_far_enough(const Body & b) override;
};

class FastPositionTracker : public PositionTracker
{
protected:
    void track_impl() override;

public:
    FastPositionTracker(const std::string & filename);
};

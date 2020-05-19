#include "Common.hpp"

#include "SimulationOutput.hpp"
#include <memory>

class OvdbConverter
{
public:
    OvdbConverter(std::shared_ptr<SimulationOutput> data);
    void Output(const std::string& path) const;

private:
    std::shared_ptr<SimulationOutput> mData;
};

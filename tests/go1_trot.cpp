#include <pinocchio/parsers/urdf.hpp>

#include "galileo/multibody/robot-spec.hpp"
#include "galileo/multibody/actuations/implementations/actuation-floating-base.hpp"
#include "galileo/multibody/states/implementations/state-multibody.hpp"

#include "galileo/predictive/phases/phase-spec.hpp"
#include "galileo/multibody/residuals/implementations/residual-frame-translation.hpp"
#include "galileo/core/activations/implementations/activation-quadratic.hpp"

#include "galileo/core/costs/implementations/cost-residual.hpp"
#include "galileo/core/costs/cost-manager.hpp"

#include "galileo/core/constraints/equality/implementations/constraint-residual.hpp"
#include "galileo/core/constraints/equality/constraint-manager.hpp"

#include "galileo/multibody/contacts/implementations/contact-3d.hpp"
#include "galileo/multibody/contacts/contact-manager.hpp"

#include "galileo/core/controls/implementations/control-param-polynomial.hpp"

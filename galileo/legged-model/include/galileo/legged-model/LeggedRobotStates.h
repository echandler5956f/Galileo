#pragma once

#include "galileo/opt/States.h"
#include "galileo/legged-model/EndEffector.h"
#include <pinocchio/autodiff/casadi.hpp>
#include "pinocchio/multibody/model.hpp"

namespace galileo
{
    namespace opt
    {
        /**
         * @brief Simple slicer class for getting state variables.
         *
         */
        class LeggedRobotStates : public States
        {
        public:
            /**
             * @brief Construct a new Legged Robot State object.
             *
             * @param nq_ The number of position variables.
             * @param nv_ The number of velocity variables.
             * @param ees The end effectors.
             */
            LeggedRobotStates(int nq_, int nv_, legged::contact::RobotEndEffectors ees);

            /**
             * @brief Get momenta: nh x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The momenta
             */
            template <class Sym>
            const Sym get_ch(const Sym &cx);

            /**
             * @brief Get momenta delta: nh x 1.
             *
             * @tparam Sym The type of the input
             * @param cdx The input
             * @return const Sym The momenta delta
             */
            template <class Sym>
            const Sym get_ch_d(const Sym &cdx);

            // /**
            //  * @brief Get momenta time derivative: nh x 1.
            //  *
            //  * @tparam Sym The type of the input
            //  * @param cx The input
            //  * @return const Sym The momenta time derivative
            //  */
            // template <class Sym>
            // const Sym get_cdh(const Sym &cx);

            // /**
            //  * @brief Get momentum time derivative delta: nh x 1.
            //  *
            //  * @tparam Sym The type of the input
            //  * @param cdx The input
            //  * @return const Sym The momentum time derivative delta.
            //  */
            // template <class Sym>
            // const Sym get_cdh_d(const Sym &cdx);

            /**/
            /**
             * @brief Get q: nq x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The position
             */
            template <class Sym>
            const Sym get_q(const Sym &cx);

            /**
             * @brief Get q delta: nv x 1.
             *
             * @tparam Sym The type of the input
             * @param cdx The input
             * @return const Sym The position delta
             */
            template <class Sym>
            const Sym get_q_d(const Sym &cdx);

            /**
             * @brief Get qj: (nq - 7) x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The joint position
             */
            template <class Sym>
            const Sym get_qj(const Sym &cx);

            // /**
            //  * @brief Get v: nv x 1.
            //  *
            //  * @tparam Sym The type of the input
            //  * @param cx The input
            //  * @return const Sym The velocity
            //  */
            // template <class Sym>
            // const Sym get_v(const Sym &cx);

            // /**
            //  * @brief Get v delta: nv x 1.
            //  *
            //  * @tparam Sym The type of the input
            //  * @param cdx The input
            //  * @return const Sym The velocity delta
            //  */
            // template <class Sym>
            // const Sym get_v_d(const Sym &cdx);

            // /**
            //  * @brief Get v_j: (nv - 6) x 1.
            //  *
            //  * @tparam Sym The type of the input
            //  * @param cx The input
            //  * @return const Sym The joint velocity
            //  */
            // template <class Sym>
            // const Sym get_vj(const Sym &cx);

            /**
             * @brief Get f: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @param ee_id The end effector id
             * @return const Sym The force
             */
            template <class Sym>
            const Sym get_f(const Sym &u, pinocchio::FrameIndex ee_id);

            /**
             * @brief Get tau: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @param ee_id The end effector id
             * @return const Sym The torque
             */
            template <class Sym>
            const Sym get_tau(const Sym &u, pinocchio::FrameIndex ee_id);

            /**
             * @brief Get the contact wrench: 6 x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @param ee_id The end effector id
             * @return const Sym The contact wrench
             */
            template <class Sym>
            const Sym get_wrench(const Sym &u, pinocchio::FrameIndex ee_id);

            /**
             * @brief Get all wrenches: nF x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @return const Sym The wrenches
             */
            template <class Sym>
            const Sym get_all_wrenches(const Sym &u);

            /**
             * @brief Get vju: nvju x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @return const Sym The joint velocity input
             */
            template <class Sym>
            const Sym get_vju(const Sym &u);

            /**
             * @brief Get the general forces: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u_general The input
             * @return const Sym The forces
             */
            template <class Sym>
            const Sym get_general_forces(const Sym &u_general);

            /**
             * @brief Get the general torques: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u_general The input
             * @return const Sym The torques
             */
            template <class Sym>
            const Sym get_general_torques(const Sym &u_general);

            /**
             * @brief Get the general joint velocities: nvju x 1.
             *
             * @tparam Sym The type of the input
             * @param u_general The input
             * @return const Sym The joint velocities
             */
            template <class Sym>
            const Sym get_general_joint_velocities(const Sym &u_general);

            // map between frame id and its index in the input vector
            std::map<pinocchio::FrameIndex, std::tuple<int, int>> frame_id_to_index_range;

            /*-----------------------------
            Sizes of certain state elements
            -----------------------------*/

            /**
             * @brief Momenta space dimension.
             *
             */
            static const int nh = 6;

            /**
             * @brief Momenta time derivative offset.
             *
             */
            static const int ndh = 6;

            /**
             * @brief Number of position coordinates for the base.
             *
             */
            static const int nqb = 7;

            /**
             * @brief Number of velocity coordinates for the base.
             *
             */
            static const int nvb = 6;

            /**
             * @brief Number of position variables.
             *
             */
            int nq = 0;

            /**
             * @brief Number of velocity variables.
             *
             */
            int nv = 0;

            /**
             * @brief Number of joint velocity inputs.
             *
             */
            int nvju = 0;

            /**
             * @brief Number of wrench variables
             *
             */
            int nF = 0;

            /**
             * @brief Summed wrenches followed by the joint velocities.
             *
             */
            int nu_general = 0;

            /*----------------------------------------
            Starting indices of certain state elements
            ----------------------------------------*/

            /**
             * @brief Starting index of the momenta.
             * 
             */
            int h_index = 0;

            // /**
            //  * @brief Starting index of the momenta time derivative.
            //  * 
            //  */
            // int dh_index = 0;

            /**
             * @brief Starting index of the position variables.
             * 
             */
            int q_index = 0;

            /**
             * @brief Starting index of the joint position variables.
             * 
             */
            int qj_index = 0;

            // /**
            //  * @brief Starting index of the velocity variables.
            //  * 
            //  */
            // int v_index = 0;

            // /**
            //  * @brief Starting index of the joint velocity inputs.
            //  * 
            //  */
            // int vj_index = 0;

            /**
             * @brief Starting index of the generalized force variables.
             * 
             */
            int general_force_index = 0;

            /**
             * @brief Starting index of the generalized torque variables.
             * 
             */
            int general_torque_index = 0;

            /**
             * @brief Starting index of the generalized joint velocity inputs.
             * 
             */
            int general_vju_index = 0;
        };

        typedef double Scalar;
        typedef casadi::SX ADScalar;

        typedef pinocchio::ModelTpl<Scalar> Model;
        typedef Model::Data Data;

        typedef pinocchio::ModelTpl<ADScalar> ADModel;
        typedef ADModel::Data ADData;

        typedef Model::ConfigVectorType ConfigVector;
        typedef Model::TangentVectorType TangentVector;

        typedef ADModel::ConfigVectorType ConfigVectorAD;
        typedef ADModel::TangentVectorType TangentVectorAD;
    }
}
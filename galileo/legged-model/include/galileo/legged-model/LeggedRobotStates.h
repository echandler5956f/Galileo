#pragma once

#include <pinocchio/autodiff/casadi.hpp>
#include "galileo/opt/States.h"

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
             * @brief Construct a new Legged Robot State object
             *
             * @param nq_nv The number of position and velocity variables.
             */
            LeggedRobotStates(int nq_, int nv_, legged::contact::RobotEndEffectors ees)
            {
                this->nq = nq_;
                this->nv = nv_;
                this->nx = this->nh + this->ndh + this->nq + this->nv;
                this->ndx = this->nh + this->ndh + 2 * this->nv;
                this->nvju = this->nv - this->nvb;
                this->nu = this->nvju;
                for (auto ee : ees)
                {
                    if (ee.second->is_6d)
                    {
                        this->frame_id_to_index_range[ee.second->frame_id] = std::make_tuple(this->nF, this->nF + 6);
                        this->nF += 6;
                    }
                    else
                    {
                        this->frame_id_to_index_range[ee.second->frame_id] = std::make_tuple(this->nF, this->nF + 3);
                        this->nF += 3;
                    }
                }
                this->nu += this->nF;
            }

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
             * @brief Get momenta: nh x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The momenta
             */
            template <class Sym>
            const Sym get_ch(const Sym &cx)
            {
                return cx(casadi::Slice(0, this->nh));
            }

            /**
             * @brief Get momenta delta: nh x 1.
             *
             * @tparam Sym The type of the input
             * @param cdx The input
             * @return const Sym The momenta delta
             */
            template <class Sym>
            const Sym get_ch_d(const Sym &cdx)
            {
                return cdx(casadi::Slice(0, this->nh));
            }

            /**
             * @brief Get momenta time derivative: nh x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The momenta time derivative
             */
            template <class Sym>
            const Sym get_cdh(const Sym &cx)
            {
                return cx(casadi::Slice(this->nh, this->nh + this->ndh));
            }

            /**
             * @brief Get momentum time derivative delta: nh x 1.
             *
             * @tparam Sym The type of the input
             * @param cdx The input
             * @return const Sym The momentum time derivative delta.
             */
            template <class Sym>
            const Sym get_cdh_d(const Sym &cdx)
            {
                return cdx(casadi::Slice(this->nh, this->nh + this->ndh));
            }

            /**/
            /**
             * @brief Get q: nq x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The position
             */
            template <class Sym>
            const Sym get_q(const Sym &cx)
            {
                return cx(casadi::Slice(this->nh + this->ndh, this->nh + this->ndh + this->nq));
            }

            /**
             * @brief Get q delta: nv x 1.
             *
             * @tparam Sym The type of the input
             * @param cdx The input
             * @return const Sym The position delta
             */
            template <class Sym>
            const Sym get_q_d(const Sym &cdx)
            {
                return cdx(casadi::Slice(this->nh + this->ndh, this->nh + this->ndh + this->nv));
            }

            /**
             * @brief Get qj: (nq - 7) x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The joint position
             */
            template <class Sym>
            const Sym get_qj(const Sym &cx)
            {
                return cx(casadi::Slice(this->nh + this->ndh + this->nqb, this->nh + this->ndh + this->nq));
            }

            /**
             * @brief Get v: nv x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The velocity
             */
            template <class Sym>
            const Sym get_v(const Sym &cx)
            {
                return cx(casadi::Slice(this->nh + this->ndh + this->nq, this->nx));
            }

            /**
             * @brief Get v delta: nv x 1.
             *
             * @tparam Sym The type of the input
             * @param cdx The input
             * @return const Sym The velocity delta
             */
            template <class Sym>
            const Sym get_v_d(const Sym &cdx)
            {
                return cdx(casadi::Slice(this->nh + this->ndh + this->nv, this->ndx));
            }

            /**
             * @brief Get v_j: (nv - 6) x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The joint velocity
             */
            template <class Sym>
            const Sym get_vj(const Sym &cx)
            {
                return cx(casadi::Slice(this->nh + this->ndh + this->nq + this->nvb, this->nx));
            }

            /**
             * @brief Get f: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @param ee_id The end effector id
             * @return const Sym The force
             */
            template <class Sym>
            const Sym get_f(const Sym &u, const int ee_id)
            {
                Sym tmp = u(casadi::Slice(this->nvju + std::get<0>(this->frame_id_to_index_range[ee_id]), this->nvju + std::get<1>(this->frame_id_to_index_range[ee_id])));
                return tmp(casadi::Slice(0, 3));
            }

            /**
             * @brief Get tau: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @param ee_id The end effector id
             * @return const Sym The torque
             */
            template <class Sym>
            const Sym get_tau(const Sym &u, const int ee_id)
            {
                Sym tmp = u(casadi::Slice(this->nvju + std::get<0>(this->frame_id_to_index_range[ee_id]), this->nvju + std::get<1>(this->frame_id_to_index_range[ee_id])));
                if (tmp.size1() < 6)
                    std::cout << "tau does not exist for ee " << ee_id << "\n";
                    return;
                else
                    return tmp(casadi::Slice(3, 6));
            }

            /**
             * @brief Get vju: nvju x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @return const Sym The joint velocity input
             */
            template <class Sym>
            const Sym get_vju(const Sym &u)
            {
                return u(casadi::Slice(this->nF, this->nu));
            }

            //map between frame id and its index in the input vector
            std::map<int, std::tuple<int, int>> frame_id_to_index_range;

            /**
             * @brief Number of position variables.
             *
             */
            int nq;

            /**
             * @brief Number of velocity variables.
             *
             */
            int nv;

            /**
             * @brief Number of joint velocity inputs.
             *
             */
            int nvju;

            /**
             * @brief Number of wrench variables
             *
             */
            int nF;
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
#pragma once

#include "galileo/legged-model/EnvironmentSurfaces.h"
#include <pinocchio/multibody/fwd.hpp>
#include <string>
#include <memory>
#include <map>

namespace galileo
{
    namespace legged
    {
        namespace contact
        {
            /**
             * @brief The types of end effectors.
             * 
             */
            enum class EE_Types
            {
                NON_PREHENSILE_3DOF,
                NON_PREHENSILE_6DOF,
                PREHENSILE_3DOF,
                PREHENSILE_6DOF
            };

            /**
             * @brief A struct for holding the end effector data of the robot.
             *
             */
            struct EndEffector
            {
                /**
                 * @brief The name of the end effector frame in pinocchio.
                 *
                 */
                std::string frame_name;

                /**
                 * @brief The id of the frame in pinocchio. Used as the key pair in global end effector maps.
                 *
                 */
                pinocchio::FrameIndex frame_id;

                /**
                 * @brief The type of the end effector.
                 *
                 */
                EE_Types ee_type;

                /**
                 * @brief The joint indices between root and this end-effector frame.
                 *
                 */
                std::vector<pinocchio::JointIndex> joint_indices;

                /**
                 * @brief Offset index in the task space vector.
                 *
                 */
                pinocchio::JointIndex offset_index;
            };

            /**
             * @brief A map-like class that uses a vector of pairs to store the data.
             *
             * Uses a vector of pairs to store the data, but provides a map-like interface ([] operator is overloaded).
             *
             */
            template <typename Key, typename Value>
            class VectorMap
            {
            public:
                /**
                 * @brief The pair type used to store the data.
                 *
                 */
                typedef std::pair<Key, Value> Pair;

                /**
                 * @brief The container type used to store the data.
                 *
                 */
                typedef std::vector<Pair> Container;

                /**
                 * @brief The iterator types.
                 *
                 */
                typedef typename Container::iterator iterator;

                /**
                 * @brief The const iterator types.
                 *
                 */
                typedef typename Container::const_iterator const_iterator;

                /**
                 * @brief Inserts a pair into the data vector.
                 *
                 * This function inserts a pair into the data vector of the EndEffector class.
                 *
                 * @param pair The pair to be inserted.
                 */
                void insert(const Pair &pair)
                {
                    data.push_back(pair);
                }

                /**
                 * @brief Overloaded [] operator to access the data.
                 *
                 * This function overloads the [] operator to access the data. If the key is not found, a new pair is inserted with the key and a default value.
                 *
                 * @param key The key to be accessed.
                 * @return Value& The value corresponding to the key.
                 */
                Value &operator[](const Key &key)
                {
                    for (Pair &pair : data)
                    {
                        if (pair.first == key)
                        {
                            return pair.second;
                        }
                    }
                    // Key not found, insert a new pair with the key and a default value
                    data.push_back(Pair(key, Value()));
                    return data.back().second;
                }

                /**
                 * @brief Accesses the value associated with the specified key.
                 *
                 * This function searches for the specified key in the data container and returns the corresponding value.
                 * If the key is found, the associated value is returned. If the key is not found, an std::out_of_range exception is thrown.
                 *
                 * @param key The key to search for.
                 * @return A reference to the value associated with the key.
                 * @throws std::out_of_range if the key is not found.
                 */
                Value &at(const Key &key)
                {
                    for (Pair &pair : data)
                    {
                        if (pair.first == key)
                        {
                            return pair.second;
                        }
                    }
                    throw std::out_of_range("Key not found");
                }

                /**
                @brief Accesses the value associated with the specified key in the data container.
                 *
                 * This function searches for the specified key in the data container and returns the associated value.
                 * If the key is not found, an std::out_of_range exception is thrown.
                 *
                 * @param key The key to search for.
                 * @return const Value& The value associated with the key.
                 * @throws std::out_of_range if the key is not found.
                 */
                const Value &at(const Key &key) const
                {
                    for (const Pair &pair : data)
                    {
                        if (pair.first == key)
                        {
                            return pair.second;
                        }
                    }
                    throw std::out_of_range("Key not found");
                }

                /**
                 * @brief Returns an iterator to the beginning of the data container.
                 *
                 * @return iterator An iterator to the beginning of the data container.
                 */
                iterator begin()
                {
                    return data.begin();
                }

                /**
                 * @brief Returns an iterator to the end of the data container.
                 *
                 * @return iterator An iterator to the end of the data container.
                 */
                iterator end()
                {
                    return data.end();
                }

                /**
                 * @brief Returns a const iterator to the beginning of the data container.
                 *
                 * @return const_iterator A const iterator to the beginning of the data container.
                 */
                const_iterator begin() const
                {
                    return data.begin();
                }

                /**
                 * @brief Returns a const iterator to the end of the data container.
                 *
                 * @return const_iterator A const iterator to the end of the data container.
                 */
                const_iterator end() const
                {
                    return data.end();
                }

                /**
                 * @brief Searches for the specified key in the data container.
                 *
                 * This function searches for the specified key in the data container and returns an iterator to the pair if the key is found.
                 * If the key is not found, an iterator to the end of the data container is returned.
                 *
                 * @param key The key to search for.
                 * @return iterator An iterator to the pair if the key is found, an iterator to the end of the data container otherwise.
                 */
                iterator find(const Key &key)
                {
                    return std::find_if(data.begin(), data.end(), [&key](const Pair &pair)
                                        { return pair.first == key; });
                }

                /**
                 * @brief Searches for the specified key in the data container.
                 *
                 * This function searches for the specified key in the data container and returns a const iterator to the pair if the key is found.
                 * If the key is not found, a const iterator to the end of the data container is returned.
                 *
                 * @param key The key to search for.
                 * @return const_iterator A const iterator to the pair if the key is found, a const iterator to the end of the data container otherwise.
                 */
                const_iterator find(const Key &key) const
                {
                    return std::find_if(data.begin(), data.end(), [&key](const Pair &pair)
                                        { return pair.first == key; });
                }

                /**
                 * @brief Returns the number of elements in the data container.
                 *
                 * @return size_t The number of elements in the data container.
                 */
                size_t size() const
                {
                    return data.size();
                }

                /**
                 * @brief Returns true if the data container is empty.
                 *
                 * @return bool True if the data container is empty, false otherwise.
                 */
                bool empty() const
                {
                    return data.empty();
                }

                /**
                 * @brief Clears the data container.
                 *
                 */
                void clear()
                {
                    data.clear();
                }

            private:
                /**
                 * @brief The underlying data container.
                 *
                 */
                Container data;
            };

            /**
             * @brief Maps the end effector index to a bool. True if the end effector is in contact.
             *
             */
            template class VectorMap<pinocchio::FrameIndex, bool>;
            typedef VectorMap<pinocchio::FrameIndex, bool> ContactCombination;

            /**
             * @brief Maps an end effector name to am End Effector ptr. We use a ptr so that there is only ever one instance of the end effector.
             *
             */
            template class VectorMap<pinocchio::FrameIndex, std::shared_ptr<EndEffector>>;
            typedef VectorMap<pinocchio::FrameIndex, std::shared_ptr<EndEffector>> RobotEndEffectors;
        }
    }
}

#ifndef FRAME_HPP
#define FRAME_HPP

#include <boost/shared_ptr.hpp>
#include "Pose.hpp"

namespace rfs{

  /** 
   * \class Frame2d
   * \brief A 2d reference frame.
   *
   * Given that the relative base frame is frame b, 
   * and the current frame is frame c,
   * this class contains the rotation information R_b_c, and translation 
   * t_b_c_c (displacement from b to c in the frame of c). If F_b points
   * to NULL, the base frame is assumed to be the inertial frame.
   */
  class Frame2d : public Pose<3, 2, 1>
  {

  public:

    /** \brief Smart pointer to Frame2d */
    typedef boost::shared_ptr<Frame2d> FramePtr;
    
    /** \brief 2d Rotation matrix */ 
    typedef Eigen::Matrix<double, 2, 2> RotMat;
    
    /** \brief Eigen vector of size 2 */
    typedef Pose<3, 2, 1>::PosVec PosVec;

    /** \brief RandomVec of size 2 */
    typedef RandomVec<2> PosRandomVec;

    /** \brief Eigen vector of size 3 */
    typedef Pose<3, 2, 1>::Vec PoseVec;
    
    /** \brief Constructor */
    Frame2d();
    
    /** 
     * \brief Constructor 
     * \param[in] T_c_b The Eigen pose vector defining this reference frame
     * \param[in] t Timestamp
     * \param[in] F_b Base frame pointer
     */
    Frame2d(PoseVec const &T_c_b, TimeStamp const &t = TimeStamp(), FramePtr F_b = FramePtr() );
    
    /** 
     * \brief Constructor 
     * \param[in] T_c_b The Pose defining this reference frame
     * \param[in] F_b Base frame pointer
     */
    Frame2d(Pose<3, 2, 1> T_c_b, FramePtr F_b = FramePtr() );
    ~Frame2d();

    /**
     * \brief Concatenate this frame F_b_c with another frame F_c_d.
     * \param[in] F_c_d other frame F_c_d
     * \return resulting frame F_b_d
     */
    Frame2d operator*(Frame2d const &F_c_d);

    /**
     * \brief Set the base frame pointer.
     * \param[in] F_b Base frame pointer
     */
    void setBaseFrame(FramePtr const &F_b);
    
    /*
     * \brief Get the rotation relative to the base frame as a matrix
     * \param[out] R_b_c Rotation matrix 
     */
    void getRotMat( RotMat &R_b_c ) const;

    /*
     * \brief Get a position in the current frame expressed in the base frame
     * \param[in,out] p_b Position of a point
     * \param[in,out] p_b_cov Covariance representing the uncertainty of p_b (optional)
     */
    void getRelToBaseFrame(PosVec &p_b, PosRandomVec::Cov *p_b_cov = NULL) const;

    /*
     * \brief Get a position in the current frame expressed in the base frame
     * \param[in] p_c Position of the point in the current frame
     * \param[out] p_b Position of the point in the base frame
     * \param[in] p_c_cov Covariance representing the uncertainty of p_c (optional)
     * \param[out] p_b_cov Covariance representing the uncertainty of p_b (optional)
     */
    void getRelToBaseFrame(PosVec const &p_c, PosVec &p_b,
			   PosRandomVec::Cov *p_c_cov = NULL, PosRandomVec::Cov *p_b_cov = NULL) const;
    /*
     * \brief Get a position in the current frame expressed in the base frame
     * \param[in,out] p_b Position of a point (which may include the covariance)
     */
    void getRelToBaseFrame(PosRandomVec &p_b) const;

    /*
     * \brief Get a position in the current frame expressed in the base frame
     * \param[in] p_c Position of the point (which may include the covariance) in the current frame 
     * \param[out] p_b Position of the point (which may include the covariance) in the base frame
     */
    void getRelToBaseFrame(PosRandomVec const &p_c, PosRandomVec &p_b) const;

  private:

    /** \brief Pointer to base reference frame */
    FramePtr F_b_;
    

  };
 
}





#endif

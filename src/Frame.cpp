#include "Frame.hpp"

namespace rfs{

  /// Implementation

  Frame2d::Frame2d(){}

  Frame2d::Frame2d(Frame2d::PoseVec const &T_c_b, TimeStamp const &t, FramePtr F_b ) :
  Pose<3, 2, 1>(T_c_b, t), F_b_(F_b){}

  Frame2d::Frame2d(Pose<3, 2, 1> T_c_b, FramePtr F_b ):
  Pose<3, 2, 1>(T_c_b), F_b_(F_b){}

  Frame2d::~Frame2d(){}

  Frame2d Frame2d::operator*(Frame2d const &F_c_d){

    double a_b_c, a_c_d;
    RotMat R_b_c, R_c_d, R_b_d;
    PosVec t_c_b_b, t_d_c_c, t_d_b_b;
    
    getRotMat(R_b_c);
    getPos(t_c_b_b);
    F_c_d.getRotMat(R_c_d);
    F_c_d.getPos(t_d_c_c);
    
    R_b_d = R_b_c * R_c_d;
    t_d_b_b = R_b_c * t_d_c_c + t_c_b_b;

    PoseVec T_b_d;
    T_b_d << t_d_b_b, atan2( R_b_d(1,0), R_b_d(0,0));
    
    TimeStamp F_b_d_time = this->getTime();
    TimeStamp F_c_d_time = F_c_d.getTime();
    if(F_b_d_time < F_c_d_time)
      F_b_d_time = F_c_d_time;
      
    return Frame2d(T_b_d, F_b_d_time, F_b_);

  }

  void Frame2d::setBaseFrame(FramePtr const &F_b){

    F_b_ = F_b;

  }

  void Frame2d::getRotMat( Frame2d::RotMat &R_b_c) const{

    double a = getRot()[0];
    R_b_c << cos(a), -sin(a), sin(a), cos(a);

  }

  void Frame2d::getRelToBaseFrame(Frame2d::PosVec &p_b, Frame2d::PosRandomVec::Cov *p_b_cov) const{

    RotMat R_b_c;
    PosVec t_c_b_b;
    getRotMat(R_b_c);
    getPos(t_c_b_b);
    
    p_b =  R_b_c * p_b + t_c_b_b;

    if(p_b_cov !=  NULL){
      *p_b_cov = R_b_c * *p_b_cov * R_b_c.transpose();
    }
  }

  void Frame2d::getRelToBaseFrame(Frame2d::PosVec const &p_c, Frame2d::PosVec &p_b,
				  Frame2d::PosRandomVec::Cov *p_c_cov, Frame2d::PosRandomVec::Cov *p_b_cov) const{

    RotMat R_b_c;
    PosVec t_c_b_b;
    getRotMat(R_b_c);
    getPos(t_c_b_b);
    
    p_b =  R_b_c * p_c + t_c_b_b;

    if(p_b_cov !=  NULL && p_b_cov !=  NULL){
      *p_b_cov = R_b_c * *p_c_cov * R_b_c.transpose();
    }
  }

  void Frame2d::getRelToBaseFrame(Frame2d::PosRandomVec &p_b) const{
    
    Frame2d::PosRandomVec::Vec p_b_tmp;
    Frame2d::PosRandomVec::Cov p_b_cov_tmp;
    p_b.get(p_b_tmp, p_b_cov_tmp);

    getRelToBaseFrame(p_b_tmp, &p_b_cov_tmp);
    p_b.set(p_b_tmp, p_b_cov_tmp);
    
  }
  
  void Frame2d::getRelToBaseFrame(Frame2d::PosRandomVec const &p_c, Frame2d::PosRandomVec &p_b) const{
    
    Frame2d::PosRandomVec::Vec p_b_tmp;
    Frame2d::PosRandomVec::Vec p_c_tmp;
    Frame2d::PosRandomVec::Cov p_b_cov_tmp;
    Frame2d::PosRandomVec::Cov p_c_cov_tmp;
    p_c.get(p_c_tmp, p_c_cov_tmp);

    getRelToBaseFrame(p_c_tmp, p_b_tmp, &p_c_cov_tmp, &p_b_cov_tmp);
    
    p_b.set(p_b_tmp, p_b_cov_tmp, p_c.getTime());

  }

}


/*
 * ---------------------------------------------------------------------------------------------------------------------
 *  Copyright 2016-2021 Doug Meyer <doug@dougdev.net>
 *
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
 *  Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 *  more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License along with this program. If not, see
 *  <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------------------------------------------------------
 */
package net.dougdev.datatype;

import static java.lang.Math.cos;
import static java.lang.Math.sin;

/**
 * Instances of this class are used to represent and manage 4x4 homogeneous transformation matrices.
 */
public final class Xform {
  private static final double[] IDENTITY = {1.0D, 0.0D, 0.0D, 0.0D,
                                            0.0D, 1.0D, 0.0D, 0.0D,
                                            0.0D, 0.0D, 1.0D, 0.0D,
                                            0.0D, 0.0D, 0.0D, 1.0D};
 
  private double[] fwd = new double[16];
  private double[] inv = new double[16];
  
  /**
   * This constructor creates a new instance and sets its value to the identity matrix.
   */
  protected Xform() { this.clear(); }
  
  /**
   * This constructor creates a new instance and sets its value to that of another specified Xform instance.
   */
  protected Xform(Xform other) { this.set(other); }
  
  /**
   * This method sets the value of the target instance to the identity matrix.
   * 
   * @return A reference to the target instance.
   */
  protected final Xform clear() {
    System.arraycopy(IDENTITY, 0, this.fwd, 0, 16);
    System.arraycopy(IDENTITY, 0, this.inv, 0, 16);
    return this;
  }
  
  /**
   * This method sets the value of the target instance to that of another specified Xform instance.
   * 
   * @return A reference to the target instance.
   */
  protected final Xform set(Xform other) {
    System.arraycopy(other.fwd, 0, this.fwd, 0, 16);
    System.arraycopy(other.inv, 0, this.inv, 0, 16);
    return this;
  }
  
  /**
   * This method performs a 3D translation operation on the target instance.
   * 
   * @param xyz An array containing the X, Y, and Z values to be used as parameters to the translation operation.
   * @return A reference to the target instance.
   */
  protected final Xform translate(double[] xyz) {
    this.fwd[12] += ((this.fwd[0] * xyz[0]) + (this.fwd[4] * xyz[1]) + (this.fwd[8] * xyz[2]));
    this.fwd[13] += ((this.fwd[1] * xyz[0]) + (this.fwd[5] * xyz[1]) + (this.fwd[9] * xyz[2]));
    this.fwd[14] += ((this.fwd[2] * xyz[0]) + (this.fwd[6] * xyz[1]) + (this.fwd[10] * xyz[2]));
    this.inv[12] -= xyz[0];
    this.inv[13] -= xyz[1];
    this.inv[14] -= xyz[2];
    return this;
  }
  
  /**
   * This method performs a 3D rotation operation on the target instance.
   * 
   * @param rpyRadians An array containing the R, P, and Y values to be used as parameters to the rotation operation.
   *                   These values are presumed to be defined in radians.
   * @return A reference to the target instance.
   */
  protected final Xform rotate(double[] rpyRadians) {
    double cosTerm, sinTerm, a, b, c, d, e, f, g, h;
    // --- Yaw angle used for rotation about Z axis ----------
    if (rpyRadians[2] != 0) {
      cosTerm = cos(rpyRadians[2]);
      sinTerm = sin(rpyRadians[2]);
      a = (this.fwd[0] * cosTerm) + (this.fwd[4] * sinTerm);
      b = (this.fwd[1] * cosTerm) + (this.fwd[5] * sinTerm);
      c = (this.fwd[2] * cosTerm) + (this.fwd[6] * sinTerm);
      d = (this.fwd[4] * cosTerm) - (this.fwd[0] * sinTerm);
      e = (this.fwd[5] * cosTerm) - (this.fwd[1] * sinTerm);
      f = (this.fwd[6] * cosTerm) - (this.fwd[2] * sinTerm);
      this.fwd[0] = a;
      this.fwd[1] = b;
      this.fwd[2] = c;
      this.fwd[4] = d;
      this.fwd[5] = e;
      this.fwd[6] = f;
      cosTerm = cos(-rpyRadians[2]);
      sinTerm = sin(-rpyRadians[2]);
      a = (this.inv[0]  * cosTerm) - (this.inv[1]  * sinTerm);
      b = (this.inv[4]  * cosTerm) - (this.inv[5]  * sinTerm);
      c = (this.inv[8]  * cosTerm) - (this.inv[9]  * sinTerm);
      d = (this.inv[12] * cosTerm) - (this.inv[13] * sinTerm);
      e = (this.inv[1]  * cosTerm) + (this.inv[0]  * sinTerm);
      f = (this.inv[5]  * cosTerm) + (this.inv[4]  * sinTerm);
      g = (this.inv[9]  * cosTerm) + (this.inv[8]  * sinTerm);
      h = (this.inv[13] * cosTerm) + (this.inv[12] * sinTerm);
      this.inv[0]  = a;
      this.inv[4]  = b;
      this.inv[8]  = c;
      this.inv[12] = d;
      this.inv[1]  = e;
      this.inv[5]  = f;
      this.inv[9]  = g;
      this.inv[13] = h;
    }
    // --- Pitch angle used for rotation about Y axis ----------
    if (rpyRadians[1] != 0) {
      cosTerm = cos(rpyRadians[1]);
      sinTerm = sin(rpyRadians[1]);
      a = (this.fwd[0]  * cosTerm) - (this.fwd[8]  * sinTerm);
      b = (this.fwd[1]  * cosTerm) - (this.fwd[9]  * sinTerm);
      c = (this.fwd[2]  * cosTerm) - (this.fwd[10] * sinTerm);
      d = (this.fwd[8]  * cosTerm) + (this.fwd[0]  * sinTerm);
      e = (this.fwd[9]  * cosTerm) + (this.fwd[1]  * sinTerm);
      f = (this.fwd[10] * cosTerm) + (this.fwd[2]  * sinTerm);
      this.fwd[0]  = a;
      this.fwd[1]  = b;
      this.fwd[2]  = c;
      this.fwd[8]  = d;
      this.fwd[9]  = e;
      this.fwd[10] = f;
      cosTerm = cos(-rpyRadians[1]);
      sinTerm = sin(-rpyRadians[1]);
      a = (this.inv[0]  * cosTerm) + (this.inv[2]  * sinTerm);
      b = (this.inv[4]  * cosTerm) + (this.inv[6]  * sinTerm);
      c = (this.inv[8]  * cosTerm) + (this.inv[10] * sinTerm);
      d = (this.inv[12] * cosTerm) + (this.inv[14] * sinTerm);
      e = (this.inv[2]  * cosTerm) - (this.inv[0]  * sinTerm);
      f = (this.inv[6]  * cosTerm) - (this.inv[4]  * sinTerm);
      g = (this.inv[10] * cosTerm) - (this.inv[8]  * sinTerm);
      h = (this.inv[14] * cosTerm) - (this.inv[12] * sinTerm);
      this.inv[0]  = a;
      this.inv[4]  = b;
      this.inv[8]  = c;
      this.inv[12] = d;
      this.inv[2]  = e;
      this.inv[6]  = f;
      this.inv[10] = g;
      this.inv[14] = h;
    }
    // --- Roll angle used for rotation about X axis -----------
    if (rpyRadians[0] != 0) {
      cosTerm = cos(rpyRadians[0]);
      sinTerm = sin(rpyRadians[0]);
      a = (this.fwd[4]  * cosTerm) + (this.fwd[8]  * sinTerm);
      b = (this.fwd[5]  * cosTerm) + (this.fwd[9]  * sinTerm);
      c = (this.fwd[6]  * cosTerm) + (this.fwd[10] * sinTerm);
      d = (this.fwd[8]  * cosTerm) - (this.fwd[4]  * sinTerm);
      e = (this.fwd[9]  * cosTerm) - (this.fwd[5]  * sinTerm);
      f = (this.fwd[10] * cosTerm) - (this.fwd[6]  * sinTerm);
      this.fwd[4]  = a;
      this.fwd[5]  = b;
      this.fwd[6]  = c;
      this.fwd[8]  = d;
      this.fwd[9]  = e;
      this.fwd[10] = f;
      cosTerm = cos(-rpyRadians[0]);
      sinTerm = sin(-rpyRadians[0]);
      a = (this.inv[1]  * cosTerm) - (this.inv[2]  * sinTerm);
      b = (this.inv[5]  * cosTerm) - (this.inv[6]  * sinTerm);
      c = (this.inv[9]  * cosTerm) - (this.inv[10] * sinTerm);
      d = (this.inv[13] * cosTerm) - (this.inv[14] * sinTerm);
      e = (this.inv[2]  * cosTerm) + (this.inv[1]  * sinTerm);
      f = (this.inv[6]  * cosTerm) + (this.inv[5]  * sinTerm);
      g = (this.inv[10] * cosTerm) + (this.inv[9]  * sinTerm);
      h = (this.inv[14] * cosTerm) + (this.inv[13] * sinTerm);
      this.inv[1]  = a;
      this.inv[5]  = b;
      this.inv[9]  = c;
      this.inv[13] = d;
      this.inv[2]  = e;
      this.inv[6]  = f;
      this.inv[10] = g;
      this.inv[14] = h;
    }
    return this;
  }
  
  /**
   * The method multiplies the values in a provided 3 element array with the forward value of the target instance. The
   * target instance is not modified by this operation.
   * 
   * @param xyz The values to be multiplied with the forward value of the target instance. This array is modified in
   *            place with the result of the operation.
   * @return A reference to the provided array argument.
   */
  protected final double[] multiply(double[] xyz) { return this.multiply(this.fwd, xyz); }
  
  /**
   * The method multiplies the values in a provided 3 element array with the inverse value of the target instance. The
   * target instance is not modified by this operation.
   * 
   * @param xyz The values to be multiplied with the inverse value of the target instance. This array is modified in
   *            place with the result of the operation.
   * @return A reference to the provided array argument.
   */
  protected final double[] multiplyInv(double[] xyz) { return this.multiply(this.inv, xyz); }
  
  /**
   * This method multiplies a specified 3 element array with a 4x4 matrix contained in a single 16 element array.
   * 
   * @param mat The 16 element array containing the 4x4 matrix. This array is not altered by this operation.
   * @param xyz The 3 element array containing the values to be multiplied. This array is modified in place with the
   *            result of the operation.
   * @return A reference to the provided 3 element array argument.
   */
  private final double[] multiply(double[] mat, double[] xyz) {
    double x, y, z, w;
    x = (mat[0] * xyz[0]) + (mat[4] * xyz[1]) + (mat[8]  * xyz[2]) + mat[12];
    y = (mat[1] * xyz[0]) + (mat[5] * xyz[1]) + (mat[9]  * xyz[2]) + mat[13];
    z = (mat[2] * xyz[0]) + (mat[6] * xyz[1]) + (mat[10] * xyz[2]) + mat[14];
    w = (mat[3] * xyz[0]) + (mat[7] * xyz[1]) + (mat[11] * xyz[2]) + mat[15];
    if (w == 1.0D) {
      xyz[0] = x;
      xyz[1] = y;
      xyz[2] = z;
    } else {
      xyz[0] = x / w;
      xyz[1] = y / w;
      xyz[2] = z / w;
    }
    return xyz;
  }
  
  /**
   * This method multiplies the target instance with the value of a second specified instance. The target instance is
   * modified in place with the result of this operation.
   * 
   * @param other The other Xform instance with which the target instance is to be multiplied.
   * @return A reference to the target instance.
   */
  protected final Xform multiply(Xform other) {
    double[] otherFwd = other.fwd;
    double[] otherInv = other.inv;
    double a, b, c, d;
    a = this.fwd[0];
    b = this.fwd[4];
    c = this.fwd[8];
    d = this.fwd[12];
    this.fwd[0]  = (a * otherFwd[0])  + (b * otherFwd[1])  + (c * otherFwd[2])  + (d * otherFwd[3]);
    this.fwd[4]  = (a * otherFwd[4])  + (b * otherFwd[5])  + (c * otherFwd[6])  + (d * otherFwd[7]);
    this.fwd[8]  = (a * otherFwd[8])  + (b * otherFwd[9])  + (c * otherFwd[10]) + (d * otherFwd[11]);
    this.fwd[12] = (a * otherFwd[12]) + (b * otherFwd[13]) + (c * otherFwd[14]) + (d * otherFwd[15]);
    a = this.fwd[1];
    b = this.fwd[5];
    c = this.fwd[9];
    d = this.fwd[13];
    this.fwd[1]  = (a * otherFwd[0])  + (b * otherFwd[1])  + (c * otherFwd[2])  + (d * otherFwd[3]);
    this.fwd[5]  = (a * otherFwd[4])  + (b * otherFwd[5])  + (c * otherFwd[6])  + (d * otherFwd[7]);
    this.fwd[9]  = (a * otherFwd[8])  + (b * otherFwd[9])  + (c * otherFwd[10]) + (d * otherFwd[11]);
    this.fwd[13] = (a * otherFwd[12]) + (b * otherFwd[13]) + (c * otherFwd[14]) + (d * otherFwd[15]);
    a = this.fwd[2];
    b = this.fwd[6];
    c = this.fwd[10];
    d = this.fwd[14];
    this.fwd[2]  = (a * otherFwd[0])  + (b * otherFwd[1])  + (c * otherFwd[2])  + (d * otherFwd[3]);
    this.fwd[6]  = (a * otherFwd[4])  + (b * otherFwd[5])  + (c * otherFwd[6])  + (d * otherFwd[7]);
    this.fwd[10] = (a * otherFwd[8])  + (b * otherFwd[9])  + (c * otherFwd[10]) + (d * otherFwd[11]);
    this.fwd[14] = (a * otherFwd[12]) + (b * otherFwd[13]) + (c * otherFwd[14]) + (d * otherFwd[15]);
    a = this.fwd[3];
    b = this.fwd[7];
    c = this.fwd[11];
    d = this.fwd[15];
    this.fwd[3]  = (a * otherFwd[0])  + (b * otherFwd[1])  + (c * otherFwd[2])  + (d * otherFwd[3]);
    this.fwd[7]  = (a * otherFwd[4])  + (b * otherFwd[5])  + (c * otherFwd[6])  + (d * otherFwd[7]);
    this.fwd[11] = (a * otherFwd[8])  + (b * otherFwd[9])  + (c * otherFwd[10]) + (d * otherFwd[11]);
    this.fwd[15] = (a * otherFwd[12]) + (b * otherFwd[13]) + (c * otherFwd[14]) + (d * otherFwd[15]);
    a = this.inv[0];
    b = this.inv[1];
    c = this.inv[2];
    d = this.inv[3];
    this.inv[0] = (a * otherInv[0]) + (b * otherInv[4]) + (c * otherInv[8])  + (d * otherInv[12]);
    this.inv[1] = (a * otherInv[1]) + (b * otherInv[5]) + (c * otherInv[9])  + (d * otherInv[13]);
    this.inv[2] = (a * otherInv[2]) + (b * otherInv[6]) + (c * otherInv[10]) + (d * otherInv[14]);
    this.inv[3] = (a * otherInv[3]) + (b * otherInv[7]) + (c * otherInv[11]) + (d * otherInv[15]);
    a = this.inv[4];
    b = this.inv[5];
    c = this.inv[6];
    d = this.inv[7];
    this.inv[4] = (a * otherInv[0]) + (b * otherInv[4]) + (c * otherInv[8])  + (d * otherInv[12]);
    this.inv[5] = (a * otherInv[1]) + (b * otherInv[5]) + (c * otherInv[9])  + (d * otherInv[13]);
    this.inv[6] = (a * otherInv[2]) + (b * otherInv[6]) + (c * otherInv[10]) + (d * otherInv[14]);
    this.inv[7] = (a * otherInv[3]) + (b * otherInv[7]) + (c * otherInv[11]) + (d * otherInv[15]);
    a = this.inv[8];
    b = this.inv[9];
    c = this.inv[10];
    d = this.inv[11];
    this.inv[8]  = (a * otherInv[0]) + (b * otherInv[4]) + (c * otherInv[8])  + (d * otherInv[12]);
    this.inv[9]  = (a * otherInv[1]) + (b * otherInv[5]) + (c * otherInv[9])  + (d * otherInv[13]);
    this.inv[10] = (a * otherInv[2]) + (b * otherInv[6]) + (c * otherInv[10]) + (d * otherInv[14]);
    this.inv[11] = (a * otherInv[3]) + (b * otherInv[7]) + (c * otherInv[11]) + (d * otherInv[15]);
    a = this.inv[12];
    b = this.inv[13];
    c = this.inv[14];
    d = this.inv[15];
    this.inv[12] = (a * otherInv[0]) + (b * otherInv[4]) + (c * otherInv[8])  + (d * otherInv[12]);
    this.inv[13] = (a * otherInv[1]) + (b * otherInv[5]) + (c * otherInv[9])  + (d * otherInv[13]);
    this.inv[14] = (a * otherInv[2]) + (b * otherInv[6]) + (c * otherInv[10]) + (d * otherInv[14]);
    this.inv[15] = (a * otherInv[3]) + (b * otherInv[7]) + (c * otherInv[11]) + (d * otherInv[15]);
    return this;
  }
  
  /**
   * This method inverts the value of the target instance. The target instance is modified in place with the result of
   * this operation.
   * 
   * @return A reference to the target instance.
   */
  protected final Xform invert() {
    double[] temp = this.fwd;
    this.fwd = this.inv;
    this.inv = temp;
    return this;
  }
  
  /**
   * This method returns the effective translation of the target instance's current value. The units of the returned
   * value is dependent on the units of values used to translate the target instance.
   * 
   * @return An array containing the effective translation of the target instance's current value.
   */
  protected final double[] getXyz() {
    double[] value = new double[3];
    value[0] = this.fwd[12];
    value[1] = this.fwd[13];
    value[2] = this.fwd[14];
    return value;
  }
  
  /**
   * This method returns the effective rotation of the target instance's current value in radians.
   * 
   * @return An array containing the effective rotation of the target instance's current value in radians.
   */
  protected final double[] getRpyRadians() {
    double[] value = new double[3];
    value[0] = Math.atan2(this.fwd[6], this.fwd[10]);
    value[1] = Math.atan2(-this.fwd[2], Math.sqrt(Math.pow(this.fwd[6], 2) + Math.pow(this.fwd[10], 2)));
    value[2] = Math.atan2(this.fwd[1], this.fwd[0]);
    return value;
  }
 
  @Override
  public final String toString() {
    return this.fwd[0] + " " + this.fwd[4] + " " + this.fwd[8]  + " " +this.fwd[12] + "\n" +
           this.fwd[1] + " " + this.fwd[5] + " " + this.fwd[9]  + " " +this.fwd[13] + "\n" +
           this.fwd[2] + " " + this.fwd[6] + " " + this.fwd[10] + " " +this.fwd[14] + "\n" +
           this.fwd[3] + " " + this.fwd[7] + " " + this.fwd[11] + " " +this.fwd[15];
  }
}

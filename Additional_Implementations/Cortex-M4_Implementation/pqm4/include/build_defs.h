/**
 *
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 2.2 (July 2025)
 *
 * Authors listed in alphabetical order:
 * 
 * @author: Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author: Marco Gianvecchio <marco.gianvecchio@mail.polimi.it>
 * @author: Patrick Karl <patrick.karl@tum.de>
 * @author: Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author: Jonas Schupp <jonas.schupp@tum.de>
 * 
 * 
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#pragma once

#define @RSDP_VARIANT@
#define CATEGORY_@category@
#define @optimiz_target@
#cmakedefine NO_TREES

/* Skip asserts so that verify() always returns a value */
#define SKIP_ASSERT 1

/******************************************************************************/
/***************************** Size optimizations *****************************/
/******************************************************************************/

/* Memory optimization flags sign operation */
#cmakedefine MEM_OPT_SIGN_INC_CMT_1
#cmakedefine MEM_OPT_SIGN_RECOMP_E_V_U
#cmakedefine MEM_OPT_SIGN_RECOMP_Y
#cmakedefine MEM_OPT_SIGN_INC_MTREE
#cmakedefine MEM_OPT_SIGN_INC_STREE
#cmakedefine MEM_OPT_SIGN_BS_SEEDS_FROM_TREE

/* Memory optimization flags verify operation */
#cmakedefine MEM_OPT_VERIFY_INC_CMT_1
#cmakedefine MEM_OPT_VERIFY_INC_Y
#cmakedefine MEM_OPT_VERIFY_INC_STREE
#cmakedefine MEM_OPT_VERIFY_INC_MTREE
#cmakedefine MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE

/* Memory optimization flags both operations */
#cmakedefine MEM_OPT_STREE_IN_CCM
#cmakedefine MEM_OPT_MTREE_IN_CCM
#cmakedefine MEM_OPT_HASH_TO_MTREE
#cmakedefine MEM_OPT_INC_SAMPLING
#cmakedefine MEM_OPT_IN_PLACE_SAMP

/* On by default:
 * - array memory reuse
 *   - sign: u y
 *   - verify: u_prime y_prime
 */

/******************************************************************************/
/***************************** Time optimizations *****************************/
/******************************************************************************/

/* Use CMSIS DSP intrinsics */
#cmakedefine TIME_OPT_SIMD

/* On by default:
 * - some arithmetic functions use lazy reduction (see AVX2-fallback impl.)
 */

/* Disabled by default:
 * - precomputed exponentiation (G^x mod P)
 */
#cmakedefine TIME_OPT_TBL_LOOKUP

/*
 * kltException.hpp
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

#ifndef KLTEXCEPTION_HPP_
#define KLTEXCEPTION_HPP_

#include <string>

namespace klt {

/** Base class for all exceptions defined in the module
 * klt.
 *
 * @ingroup klt
 */
class KltException{

public:

	/** This enumeration defines exceptions id for the module
	 * klt.
	 */
	enum ExceptionId {
		KLT_ERROR
		//        MY_ERROR /**< my error */
	};

	/** Constructor. You should not use this constructor directly,
	 * prefer macros jfrThrowEx or jfrCreateEx which fill for you
	 * parameters \c file_ and \c line_.
	 *
	 * @param id_ exception id
	 * @param message_ message used for debug
	 * @param file_ where the exception was thrown
	 * @param line_ where the exception was thrown
	 */
	KltException(ExceptionId id_,
			const std::string& message_,
			const std::string& file_, int line_) throw();

	virtual ~KltException() throw();

	ExceptionId getExceptionId() const throw();

protected:

	ExceptionId id;

	static std::string exceptionIdToString(ExceptionId id_) throw();

}; // class KltException

} // namespace klt

#endif /* KLTEXCEPTION_HPP_ */

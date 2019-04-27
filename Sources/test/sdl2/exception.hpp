#ifndef SDL2_EXCEPTION_HPP
#define SDL2_EXCEPTION_HPP

#include <exception>
#include <string>

namespace SDL2 {
    class Exception : public std::exception {
        public:
            explicit Exception(const std::string& msg = "", bool getSDlError = true) noexcept;

            Exception(const Exception& e) noexcept {
                m_What = e.m_What;
            };

            ~Exception() override = default;

            const char* what() const noexcept override;

        protected:
            std::string m_What;
    };

    class VKException : public Exception {
        public:
            explicit VKException(const std::string& msg = "", bool getSDlError = false) noexcept;

            ~VKException() override = default;
    };
}

#endif

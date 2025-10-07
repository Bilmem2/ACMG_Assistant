"""
API Error Handler Module
=========================

Centralized error handling, logging, and retry logic for all API calls.
Provides consistent error messages, fallback strategies, and monitoring.
"""

import logging
import time
import json
from typing import Dict, Any, Optional, Callable, Tuple
from datetime import datetime
from enum import Enum
from colorama import Fore, Style, init
import requests

# Initialize colorama
init()


class ErrorSeverity(Enum):
    """Error severity levels for API failures."""
    CRITICAL = "CRITICAL"  # API down, no fallback available
    HIGH = "HIGH"          # API error but fallback exists
    MEDIUM = "MEDIUM"      # Temporary failure, retry possible
    LOW = "LOW"            # Minor issue, degraded functionality
    INFO = "INFO"          # Informational, not an error


class APIErrorType(Enum):
    """Types of API errors."""
    NETWORK_ERROR = "Network connection failed"
    TIMEOUT_ERROR = "Request timed out"
    RATE_LIMIT_ERROR = "Rate limit exceeded"
    AUTH_ERROR = "Authentication failed"
    NOT_FOUND_ERROR = "Resource not found"
    SERVER_ERROR = "Server error (5xx)"
    CLIENT_ERROR = "Client error (4xx)"
    PARSE_ERROR = "Failed to parse response"
    VALIDATION_ERROR = "Response validation failed"
    UNKNOWN_ERROR = "Unknown error"


class APIErrorHandler:
    """
    Centralized error handling for all API calls.
    
    Features:
    - Automatic retry with exponential backoff
    - Error logging with context
    - Fallback strategy execution
    - Error statistics tracking
    - User-friendly error messages
    """
    
    def __init__(self, log_file: str = "api_errors.log", enable_console: bool = True):
        """
        Initialize error handler.
        
        Args:
            log_file (str): Path to error log file
            enable_console (bool): Whether to print errors to console
        """
        self.log_file = log_file
        self.enable_console = enable_console
        self.error_stats = {
            'total_errors': 0,
            'errors_by_api': {},
            'errors_by_type': {},
            'last_error_time': None
        }
        
        # Setup logging
        self._setup_logging()
    
    def _setup_logging(self):
        """Configure logging for API errors."""
        self.logger = logging.getLogger('APIErrorHandler')
        self.logger.setLevel(logging.DEBUG)
        
        # File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        
        # Formatter
        formatter = logging.Formatter(
            '%(asctime)s | %(levelname)-8s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        
        self.logger.addHandler(file_handler)
    
    def handle_api_call(
        self,
        api_name: str,
        func: Callable,
        max_retries: int = 3,
        retry_delay: float = 1.0,
        exponential_backoff: bool = True,
        fallback_func: Optional[Callable] = None,
        fallback_return: Optional[Any] = None,
        **kwargs
    ) -> Tuple[Any, Optional[str]]:
        """
        Execute API call with error handling and retry logic.
        
        Args:
            api_name (str): Name of API being called
            func (Callable): Function to execute
            max_retries (int): Maximum retry attempts
            retry_delay (float): Initial delay between retries (seconds)
            exponential_backoff (bool): Use exponential backoff for retries
            fallback_func (Callable, optional): Fallback function if all retries fail
            fallback_return (Any, optional): Default return value if all fails
            **kwargs: Arguments to pass to func
        
        Returns:
            Tuple[Any, Optional[str]]: (result, error_message)
        """
        last_error = None
        error_type = None
        
        for attempt in range(max_retries + 1):
            try:
                # Execute the API call
                result = func(**kwargs)
                
                # Success - log if this was a retry
                if attempt > 0:
                    self._log_success_after_retry(api_name, attempt)
                
                return result, None
            
            except requests.Timeout as e:
                error_type = APIErrorType.TIMEOUT_ERROR
                last_error = str(e)
                severity = ErrorSeverity.MEDIUM
            
            except requests.ConnectionError as e:
                error_type = APIErrorType.NETWORK_ERROR
                last_error = str(e)
                severity = ErrorSeverity.HIGH
            
            except requests.HTTPError as e:
                if e.response.status_code == 429:
                    error_type = APIErrorType.RATE_LIMIT_ERROR
                    severity = ErrorSeverity.MEDIUM
                elif 500 <= e.response.status_code < 600:
                    error_type = APIErrorType.SERVER_ERROR
                    severity = ErrorSeverity.HIGH
                elif 400 <= e.response.status_code < 500:
                    error_type = APIErrorType.CLIENT_ERROR
                    severity = ErrorSeverity.LOW
                else:
                    error_type = APIErrorType.UNKNOWN_ERROR
                    severity = ErrorSeverity.MEDIUM
                last_error = f"HTTP {e.response.status_code}: {str(e)}"
            
            except json.JSONDecodeError as e:
                error_type = APIErrorType.PARSE_ERROR
                last_error = f"JSON parse error: {str(e)}"
                severity = ErrorSeverity.LOW
            
            except Exception as e:
                error_type = APIErrorType.UNKNOWN_ERROR
                last_error = str(e)
                severity = ErrorSeverity.MEDIUM
            
            # Log the error
            self._log_error(
                api_name=api_name,
                error_type=error_type,
                error_message=last_error,
                severity=severity,
                attempt=attempt + 1,
                max_retries=max_retries,
                kwargs=kwargs
            )
            
            # Update statistics
            self._update_error_stats(api_name, error_type)
            
            # Check if we should retry
            if attempt < max_retries:
                # Calculate delay with exponential backoff
                if exponential_backoff:
                    delay = retry_delay * (2 ** attempt)
                else:
                    delay = retry_delay
                
                # Log retry attempt
                self._log_retry(api_name, attempt + 1, max_retries, delay)
                time.sleep(delay)
            else:
                # All retries exhausted
                self._log_all_retries_failed(api_name, max_retries, last_error)
        
        # All retries failed - try fallback
        if fallback_func is not None:
            try:
                self._log_fallback_attempt(api_name)
                result = fallback_func(**kwargs)
                return result, f"Using fallback (API error: {last_error})"
            except Exception as e:
                self._log_fallback_failed(api_name, str(e))
        
        # Return default fallback value
        if fallback_return is not None:
            return fallback_return, f"API failed, using default value (error: {last_error})"
        
        # Complete failure
        if error_type is None:
            error_type = APIErrorType.UNKNOWN_ERROR
        if last_error is None:
            last_error = "Unknown error occurred"
        
        error_message = self._format_final_error_message(api_name, error_type, last_error)
        return None, error_message
    
    def _log_error(
        self,
        api_name: str,
        error_type: APIErrorType,
        error_message: str,
        severity: ErrorSeverity,
        attempt: int,
        max_retries: int,
        kwargs: Dict[str, Any]
    ):
        """Log error with full context."""
        log_message = (
            f"API: {api_name} | "
            f"Type: {error_type.name} | "
            f"Severity: {severity.value} | "
            f"Attempt: {attempt}/{max_retries + 1} | "
            f"Error: {error_message} | "
            f"Args: {self._sanitize_kwargs(kwargs)}"
        )
        
        # Log to file
        if severity in [ErrorSeverity.CRITICAL, ErrorSeverity.HIGH]:
            self.logger.error(log_message)
        elif severity == ErrorSeverity.MEDIUM:
            self.logger.warning(log_message)
        else:
            self.logger.info(log_message)
        
        # Console output
        if self.enable_console:
            color = self._get_severity_color(severity)
            icon = self._get_severity_icon(severity)
            print(f"{color}{icon} {api_name} Error ({error_type.name}): {error_message}{Style.RESET_ALL}")
    
    def _log_retry(self, api_name: str, attempt: int, max_retries: int, delay: float):
        """Log retry attempt."""
        message = f"API: {api_name} | Retrying in {delay:.1f}s (attempt {attempt}/{max_retries})"
        self.logger.info(message)
        
        if self.enable_console:
            print(f"{Fore.YELLOW}🔄 Retrying {api_name} in {delay:.1f}s...{Style.RESET_ALL}")
    
    def _log_success_after_retry(self, api_name: str, attempts: int):
        """Log successful call after retry."""
        message = f"API: {api_name} | Success after {attempts} retries"
        self.logger.info(message)
        
        if self.enable_console:
            print(f"{Fore.GREEN}✅ {api_name} succeeded after {attempts} retries{Style.RESET_ALL}")
    
    def _log_all_retries_failed(self, api_name: str, max_retries: int, error: str):
        """Log complete failure after all retries."""
        message = f"API: {api_name} | All {max_retries} retries failed | Last error: {error}"
        self.logger.error(message)
        
        if self.enable_console:
            print(f"{Fore.RED}❌ {api_name} failed after {max_retries} retries{Style.RESET_ALL}")
    
    def _log_fallback_attempt(self, api_name: str):
        """Log fallback strategy attempt."""
        message = f"API: {api_name} | Attempting fallback strategy"
        self.logger.info(message)
        
        if self.enable_console:
            print(f"{Fore.CYAN}🔀 Using fallback for {api_name}...{Style.RESET_ALL}")
    
    def _log_fallback_failed(self, api_name: str, error: str):
        """Log fallback failure."""
        message = f"API: {api_name} | Fallback failed | Error: {error}"
        self.logger.error(message)
        
        if self.enable_console:
            print(f"{Fore.RED}❌ Fallback for {api_name} also failed: {error}{Style.RESET_ALL}")
    
    def _update_error_stats(self, api_name: str, error_type: APIErrorType):
        """Update error statistics."""
        self.error_stats['total_errors'] += 1
        self.error_stats['last_error_time'] = datetime.now().isoformat()
        
        # Count by API
        if api_name not in self.error_stats['errors_by_api']:
            self.error_stats['errors_by_api'][api_name] = 0
        self.error_stats['errors_by_api'][api_name] += 1
        
        # Count by type
        error_type_name = error_type.name
        if error_type_name not in self.error_stats['errors_by_type']:
            self.error_stats['errors_by_type'][error_type_name] = 0
        self.error_stats['errors_by_type'][error_type_name] += 1
    
    def get_error_statistics(self) -> Dict[str, Any]:
        """
        Get error statistics.
        
        Returns:
            Dict with error counts by API and type
        """
        return self.error_stats.copy()
    
    def print_error_summary(self):
        """Print error summary to console."""
        stats = self.error_stats
        
        print(f"\n{Fore.CYAN}{'='*80}{Style.RESET_ALL}")
        print(f"{Fore.CYAN}API ERROR SUMMARY{Style.RESET_ALL}")
        print(f"{Fore.CYAN}{'='*80}{Style.RESET_ALL}\n")
        
        print(f"Total Errors: {stats['total_errors']}")
        print(f"Last Error: {stats['last_error_time'] or 'None'}\n")
        
        if stats['errors_by_api']:
            print(f"{Fore.YELLOW}Errors by API:{Style.RESET_ALL}")
            for api, count in sorted(stats['errors_by_api'].items(), key=lambda x: x[1], reverse=True):
                print(f"  {api}: {count}")
            print()
        
        if stats['errors_by_type']:
            print(f"{Fore.YELLOW}Errors by Type:{Style.RESET_ALL}")
            for error_type, count in sorted(stats['errors_by_type'].items(), key=lambda x: x[1], reverse=True):
                print(f"  {error_type}: {count}")
            print()
    
    @staticmethod
    def _get_severity_color(severity: ErrorSeverity) -> str:
        """Get console color for severity level."""
        colors = {
            ErrorSeverity.CRITICAL: Fore.RED + Style.BRIGHT,
            ErrorSeverity.HIGH: Fore.RED,
            ErrorSeverity.MEDIUM: Fore.YELLOW,
            ErrorSeverity.LOW: Fore.CYAN,
            ErrorSeverity.INFO: Fore.WHITE
        }
        return colors.get(severity, Fore.WHITE)
    
    @staticmethod
    def _get_severity_icon(severity: ErrorSeverity) -> str:
        """Get icon for severity level."""
        icons = {
            ErrorSeverity.CRITICAL: "🚨",
            ErrorSeverity.HIGH: "❌",
            ErrorSeverity.MEDIUM: "⚠️",
            ErrorSeverity.LOW: "ℹ️",
            ErrorSeverity.INFO: "📝"
        }
        return icons.get(severity, "❓")
    
    @staticmethod
    def _sanitize_kwargs(kwargs: Dict[str, Any]) -> str:
        """Sanitize kwargs for logging (remove sensitive data)."""
        # Remove potentially large or sensitive fields
        sanitized = {}
        for key, value in kwargs.items():
            if key in ['password', 'token', 'api_key', 'secret']:
                sanitized[key] = '***REDACTED***'
            elif isinstance(value, (str, int, float, bool)):
                sanitized[key] = value
            elif value is None:
                sanitized[key] = None
            else:
                sanitized[key] = f"<{type(value).__name__}>"
        
        return str(sanitized)
    
    @staticmethod
    def _format_final_error_message(
        api_name: str,
        error_type: APIErrorType,
        error_message: str
    ) -> str:
        """Format user-friendly final error message."""
        return (
            f"{api_name} API failed: {error_type.value}. "
            f"Details: {error_message}. "
            f"Please check your internet connection or try again later."
        )


# Singleton instance
_error_handler_instance = None


def get_error_handler() -> APIErrorHandler:
    """
    Get singleton error handler instance.
    
    Returns:
        APIErrorHandler: Global error handler instance
    """
    global _error_handler_instance
    if _error_handler_instance is None:
        _error_handler_instance = APIErrorHandler()
    return _error_handler_instance

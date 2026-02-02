

/// Error type for `thor`  
/// 
/// WIP: perhaps this should be multiple errors?
#[derive(Debug)]
pub enum ThorError {
    FileOpenError, 
    FloatParseError, 
    ColNumberChange
}
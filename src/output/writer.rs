use crate::myth::Myth;

use crate::error::Result;

pub trait MythWriter {
    fn write_myth(&mut self, myth: Myth) -> Result<()>;
    fn end_batch(&mut self) -> Result<()>;
    fn close(self) -> Result<()>;
}
